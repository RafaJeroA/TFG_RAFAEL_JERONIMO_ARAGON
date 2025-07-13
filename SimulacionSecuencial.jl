module SimulacionSecuencial
using ..Datos, ..Agentes, ..Tecnologias, ..PlanificacionMes, ..Resultados, ..OptimizacionSistema, ..Escenarios
using ..PoliticaVerde
using ..Resultados: append_snapshot!, exportar_escenario, graficar_snapshots, graficar_traza, df_snapshots
using Graphs
using CSV, DataFrames  
import Dates: Date, year, month, Month, Day, day
using MathOptInterface, Distributions
using ..InversionesPendientes
const MOI = MathOptInterface
export ejecutar_mercado_real  

using Base.Threads

include("EquilibrioMultivector.jl")
using .EquilibrioMultivector


const FULL_AUTOCONSUMO = Dict{Date,Dict{String,Dict{Int,Float64}}}()
const produccion_mensual = Dict{Int,Dict{Int,Dict{String,Float64}}}()
const deficits_mensual   = Dict{Int,Dict{Int,Dict{String,Float64}}}()

const LAST_LCP_SOLUTION = Dict{Date,Vector{Float64}}()

const ACTIVACION_PREDEF_MILP = false # Activar para que utilice el MILP como prioritario y PlanificaciónMes como falloback. Como el MILP tarda mucho y da algunos problemas lo he desactivado
#–– Registro de inversiones pendientes (mes_listo, tecnología, mw) por agente




"""
    ejecutar_mercado_real(escenario, agentes_iniciales, datos_bundle, tecnologias)

Ejecuta una simulación secuencial del mercado energético con decisiones mensuales de inversión
y despacho diario para todos los vectores energéticos.

# Argumentos
- `escenario::Escenario`: Configuración del escenario (nombre, políticas, etc.)
- `agentes_iniciales::Dict{Int,EnergyAgent}`: Diccionario de agentes energéticos
- `datos_bundle::NamedTuple`: Bundle con datos del sistema (demanda, perfiles, políticas, etc.)
- `tecnologias::Dict{String,Tecnologia}`: Diccionario de tecnologías disponibles

# Retorna
NamedTuple con resultados de la simulación incluyendo capacidades, producción, finanzas y precios.
"""
function ejecutar_mercado_real(escenario, agentes_iniciales, datos_bundle, tecnologias)
    # ═══════════════════════════════════════════════════════════════════════
    # VALIDACIÓN DE INPUTS
    # ═══════════════════════════════════════════════════════════════════════
    validar_datos_entrada(datos_bundle, tecnologias)
    empty!(FULL_AUTOCONSUMO)
    empty!(produccion_mensual)
    empty!(deficits_mensual)
    InversionesPendientes.limpiar_pendientes!()

    finanzas_anuales = Dict{Int, Vector{NamedTuple}}()
    println("🚀 Iniciando simulación de mercado real para $(escenario.nombre)")
    
    # Copiar agentes para no modificar originales
    agentes = copy(agentes_iniciales)
    
    # Plantilla base para las decisiones de inversión
    capex_base = OptimizacionSistema.build_model_capex_base(
        agentes,
        tecnologias,
        horizon_months = 60
    )
    
    # Extraer años desde datos_bundle.demanda
    anios = sort(unique(Int.(datos_bundle.demanda.anio)))
    println("📅 Simulando años: $(first(anios))-$(last(anios))")
    
    # Histórico de precios para decisiones de inversión
    precios_historicos = Dict{ Date, Dict{ Tuple{Int,String}, Float64 } }()

    # Resultados agregados
    capacidades_finales = Dict{Int,Dict{String,Float64}}()
    produccion_anual = Dict{Int,Dict{String,Float64}}()
    


    # ═══════════════════════════════════════════════════════════════════════
    # BUCLE PRINCIPAL POR AÑOS
    # ═══════════════════════════════════════════════════════════════════════
    for anio in anios
        for ag in values(agentes)
            if hasproperty(ag, :produccion_anual)
                ag.produccion_anual = Dict{String,Float64}()
            end
        end
        println("\n📆 Procesando año: $anio")

        # --- APLICAR RETIRO PROGRAMADO DE CAPACIDAD ---
        if hasproperty(escenario, :retiro_programado) && !isempty(escenario.retiro_programado)
            for (tec_a_retirar, retiros_anio) in escenario.retiro_programado
                if haskey(retiros_anio, anio)
                    mw_total = retiros_anio[anio]
                    println(" Retiro programado para '$tec_a_retirar' en $anio: $(mw_total) MW")

                    # 1) Capacidad total existente
                    cap_total = sum(get(ag.capacidades, tec_a_retirar, 0.0) for ag in values(agentes))
                    if cap_total > 0
                        # 2) Retiro proporcional por agente
                        for ag in values(agentes)
                            cap_ag = get(ag.capacidades, tec_a_retirar, 0.0)
                            if cap_ag > 0
                                propor = cap_ag / cap_total
                                mw_ag = mw_total * propor
                                ag.capacidades[tec_a_retirar] = max(0.0, cap_ag - mw_ag)
                                println("  ⚠️   - $(ag.name): retira $(round(mw_ag,digits=2)) MW de producción vía $tec_a_retirar; queda $(round(ag.capacidades[tec_a_retirar],digits=2)) MW")

                                # 3) Coste de desmantelamiento (€100k/MW)
                                coste = mw_ag * 100_000
                                ag.cash -= coste
                                ag.cost += coste
                            end
                        end
                    else
                        println("  ⚠️   - No hay capacidad de '$tec_a_retirar' para retirar.")
                    end
                end
            end
        end
        # --- FIN RETIRO PROGRAMADO ---
        produccion_anual[anio] = Dict{String, Float64}()
        produccion_mensual[anio] = Dict{Int,Dict{String,Float64}}()
        deficits_mensual[anio]   = Dict{Int,Dict{String,Float64}}()
        subv_unit_mes = Dict{Int,Dict{String,Float64}}()
        empty!(subv_unit_mes)

        # Centralizar extracción de políticas para el año
        politicas_anio = extraer_politicas_anio(datos_bundle.politicas, anio)
        precio_co2_anio = politicas_anio.precio_co2
        subv_max_anual = politicas_anio.subv_max_anual  # subvención máxima anual, €
        subv_max_restante_anual =  subv_max_anual        # saldo del año
        
        deficits_reales_por_vector = Dict{String,Float64}()  # Deficit real del mes en curso, clave = vector (REVISAR, NO DEBERÍA METERLO DENTRO DEL BUCLE MESES??)

        for mes in 1:12
            println("  📅 Procesando mes: $mes")
            
            # Acumulará el déficit MWh del mes en curso, clave = vector


            # ── Preparar demanda mensual agregada ──
            fechas_agregacion = Datos.fechas_del_mes(anio, mes)
            demanda_mes_df = agregar_demanda_mensual(datos_bundle, fechas_agregacion)
            
            # ── Generar perfil de viento aleatorio ──
            exp_eolica = generar_perfil_eolico_aleatorio(size(demanda_mes_df, 1))
            
            # Aplicar nuevas capacidades tras completar el delay
            for ag in values(agentes)
                # "aplicadas" ya devuelve [(tecnologia, mw), …] si inv[1] <= mes
                aplicadas = InversionesPendientes.aplicar_pendientes!(ag.id, anio*12 + mes)
                for (tec, mw) in aplicadas
                    ag.capacidades[tec] = get(ag.capacidades, tec, 0.0) + mw
                    # ── NUEVO ── asegurar que el agente conoce los costes
                    #             variables y las emisiones de la tecnología
                    if !haskey(ag.var_costs, tec)
                        ag.var_costs[tec] = tecnologias[tec].costo_om_variable
                    end
                    if !haskey(ag.emissions, tec)
                        ag.emissions[tec] = tecnologias[tec].emisiones_unitarias
                    end                    
                    # 2. Actualizar la capacidad del sistema
                    periodo = (anio, mes)
                    cap_dict = get!(OptimizacionSistema.cap_existente, tec,
                                    Dict{OptimizacionSistema.Periodo,Float64}())
                    cap_dict[periodo] = get(cap_dict, periodo, 0.0) + mw
                    # Registrar el gasto en capex justo al activarse
                    capex = mw * tecnologias[tec].costo_inicial
                    Agentes.registrar_inversion!(ag, capex)
                    println("    ✅ Activada capacidad: $(ag.name) - $tec: $mw mw (capex = $capex €)")
                end
                if !isempty(aplicadas)
                    println("    📈 Capacidades actualizadas de $(ag.name):")
                    for (tec, cap) in ag.capacidades
                        if cap > 1e-6
                            println("      - $tec: $(round(cap, digits=1)) MW")
                        end
                    end
                end
            end
            # Guardar déficits del día para decisiones de inversión
            # El déficit real del mes anterior llega calculado desde el mercado diario
            deficits_acumulados = copy(deficits_reales_por_vector)

            get!(deficits_mensual[anio], mes, Dict{String,Float64}())
            # ─ fin de mes ─
            deficits_mensual[anio][mes] =
                    copy(deficits_reales_por_vector)   # guarda lo que ocurrió
            empty!(deficits_reales_por_vector)         # resetea para el mes siguiente



            # Reiniciar acumulador para el mes que comienza
            empty!(deficits_reales_por_vector)

            #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

            # ── DÍA 1 DEL MES: DECISIONES DE INVERSIÓN ──
            if (mes - 1) % 3 == 0  # enero, abril, julio, octubre …
                # 👉 En el primer mes todavía NO hay precios históricos → nada de inversiones
                if isempty(precios_historicos)
                    println("    ⏸️ Inversiones pospuestas hasta contar con precios reales.")
                else
                    # Resetear contadores mensuales de todos los agentes
                    for ag in values(agentes)
                        ag.profit_mes = 0.0
                        ag.emissions_mes = 0.0
                    end

                    # ──────────────────────────────────────────────────────
                    # Media de precios (€/MWh) por vector a partir de los
                    # últimos 30 días y réplica a cada tecnología
                    # ──────────────────────────────────────────────────────
                    precio_medio_vector = Dict{String,Float64}()
                    if !isempty(precios_historicos)
                        acum = Dict{String,Tuple{Float64,Int}}()
                        for f in sort(collect(keys(precios_historicos)), rev=true)[1:min(30, length(precios_historicos))]
                            for ((h, vec), p) in precios_historicos[f]
                                s, n = get(acum, vec, (0.0, 0))
                                acum[vec] = (s+p, n+1)
                            end
                        end
                        for (vec, (s, n)) in acum
                            precio_medio_vector[vec] = s / n
                        end
                        # Copiar la media del vector a cada tecnología
                        for (tec_name, tec) in tecnologias
                            precio_medio_vector[tec_name] =
                                get(precio_medio_vector, tec.vector, 0.0)
                        end
                    end
                    println("    💰 Evaluando inversiones estratégicas para el trimestre…")
                    invertir_mensualmente_con_milp!(
                        agentes,
                        tecnologias,
                        datos_bundle,
                        escenario,
                        anio,
                        mes,
                        precios_historicos,
                        deficits_acumulados,
                        precio_co2_anio,
                        precio_medio_vector
                    )
                end
            end
            fecha_mes = Date(anio, mes, 1)
            #precios_mes = get(precios_historicos, fecha_mes, Dict{ Tuple{Int,String,String,Int}, Float64 }())
            precios_mes = precios_historicos   # sigue guardando precios para reporting

            # ------------------------- FIN de la zona modificada -------------------------
            



             
            # NO ME ACUERDO DE CÓMO SE LLAMABA EL PERFIL DE GENERACIÓN, POR LO QUE LO VUELVO A CONVOCAR
            # Lo convoco aquí para que no se recalcule cada día
            perfiles_gen_file = joinpath("Datos","perfiles_horarios.csv")
            df_hor_generacion = DataFrame(tecnologia=String[], estacion=String[], hora=Int[], factor=Float64[])
            try
                if !isfile(perfiles_gen_file) || filesize(perfiles_gen_file) == 0
                    @warn "'perfiles_horarios.csv' inexistente o vacío. Productores usarán factor 1.0."
                else
                    df_hor_generacion = CSV.read(perfiles_gen_file, DataFrame)
                    normalizar_cabeceras!(df_hor_generacion)
                    @info "Perfiles horarios de GENERACIÓN cargados desde $(basename(perfiles_gen_file))"
                end
            catch e
                @warn "No se pudo cargar perfiles_horarios.csv: $e. Productores usarán factor 1.0." 
            end



            # ─────────────────────────────────────────────────────────────
            # Capacidad horaria (MW) por (agente, tecnología, hora)
            # ─────────────────────────────────────────────────────────────
            capacidad_disponible = Dict{Tuple{Int,Symbol,Int},Float64}()

            for (ag_id, ag) in agentes
                for (tec, mw) in ag.capacidades
                    perfil = Dict{Int,Float64}(h => 1.0 for h in 0:23)

                    if !isempty(df_hor_generacion)
                        filas = filter(r -> r.tecnologia == tec, df_hor_generacion)
                        for row in eachrow(filas)
                            if !ismissing(row.factor) && isfinite(row.factor)
                                perfil[row.hora] = Float64(row.factor)
                            end
                        end
                    end
                    for h in 0:23
                        capacidad_disponible[(ag_id, Symbol(tec), h)] = mw * perfil[h] #* perf_al[tec][h]
                    end
                end
            end
            # cierre recalculo perf diario
            


            # ── DÍAS DEL MES: SUBASTAS DIARIAS ──
            fechas_subastas = Datos.fechas_del_mes(anio, mes)
            for fecha in fechas_subastas
                # Procesar mercado diario integrado
                subv = get(subv_unit_mes, mes, Dict{String,Float64}())
                resp = procesar_mercado_diario_integrado!(
                    fecha, escenario, agentes, datos_bundle,
                    tecnologias, precio_co2_anio, produccion_anual,
                    produccion_mensual,
                    subv_unit_mes, mes, capacidad_disponible
                )
                # Sumar los déficits físicos de hoy al acumulador del mes
                for (v, mwh) in resp.deficits
                    deficits_reales_por_vector[v] =
                        get(deficits_reales_por_vector, v, 0.0) + mwh
                end
                # Guardar precios históricos
                precios_historicos[fecha] = resp.precios
            end
            
            # ── FIN DE MES: CONTABILIDAD ──
            for ag in values(agentes)
                Agentes.actualizar_finanzas_mes!(ag, anio, mes)
            end
            Agentes.pagar_cuotas_mensuales!(agentes)

            # ── Política verde mensual ────────────────────────────────────────────────
            benef_mes = Dict(ag.id => ag.profit_mes for ag in values(agentes))

            gen_verde_mes = Dict{String,Float64}()
            for (tec, mwh) in produccion_anual[anio]          # ajustar si dicc. difiere
                if tecnologias[tec].verde
                    gen_verde_mes[tec] = get(gen_verde_mes, tec, 0.0) + mwh
                end
            end

            emisiones_mes = sum(ag.emissions_mes for ag in values(agentes))

            subv_unit = PoliticaVerde.fijar_politica!(
                anio, mes,
                emisiones_mes,
                benef_mes,
                precio_co2_anio,
                subv_max_restante_anual,
                gen_verde_mes
            )

            subv_unit_mes[mes] = subv_unit
            # Actualiza el saldo anual (subv_unit es negativo → ingreso para el agente)
            #–– Calcular subvención total del mes y descontarla del presupuesto anual ––
            total_subv_mes = sum(subv_unit[tec] * get(gen_verde_mes, tec, 0.0) for tec in keys(subv_unit); init=0.0)
            subv_max_restante_anual = max(subv_max_restante_anual - total_subv_mes, 0.0)
            # Si ya no queda presupuesto, en meses siguientes subv_unit será cero automáticamente            
            # ───────────────────────────────────────────────────────────────────────────


        end
        
        # ── FIN DE AÑO: AGREGAR RESULTADOS ──
        finanzas_anuales[anio] = construir_snapshot_finanzas(agentes, anio)
        capacidades_finales[anio] = agregar_capacidades_sistema(agentes, tecnologias)
        
        # Agregar producción anual de agentes
        for ag in values(agentes)
            if hasproperty(ag, :produccion_anual)
                for (tec, mwh) in ag.produccion_anual
                    println("      📊 Agente $(ag.name) produjo $mwh MWh con $tec")  # Debug
                    produccion_anual[anio][tec] = get(produccion_anual[anio], tec, 0.0) + mwh
                end
            end
        end

        # ─── Resumen anual de producción ───────────────────────────
        if haskey(produccion_anual, anio)
            println("\n📊  RESUMEN PRODUCCIÓN AÑO $anio")
            println("="^60)
            total = 0.0
            for (tec, mwh) in sort(collect(produccion_anual[anio]); by = x->x[1])
                println(rpad(tec, 30) *
                        lpad(string(round(mwh/1e6, digits = 2)), 10) * "  TWh")
                total += mwh
            end
            println("-"^60)
            println(rpad("TOTAL", 30) *
                    lpad(string(round(total/1e6, digits = 2)), 10) * "  TWh")
            println("="^60)
        end

        # Resetear contadores anuales
        for ag in values(agentes)
            ag.revenues = 0.0
            ag.cost = 0.0
            ag.profit = 0.0
        end
    end
    
    # ═══════════════════════════════════════════════════════════════════════
    # CONSTRUIR RESULTADOS FINALES
    # ═══════════════════════════════════════════════════════════════════════
    #=return construir_resultados_simulacion(
        escenario, agentes, capacidades_finales,
        produccion_anual, finanzas_anuales, precios_historicos,
        produccion_mensual, deficits_mensual, datos_bundle
    )
end=#

    return construir_resultados_simulacion(
            escenario,
            agentes,
            capacidades_finales,
            produccion_anual,
            finanzas_anuales,
            precios_historicos,
            produccion_mensual,
            deficits_mensual,
            datos_bundle
        )
end

"""
    validar_datos_entrada(datos_bundle::NamedTuple, tecnologias::Dict)

Valida que los datos de entrada contengan todas las columnas esperadas
y cubran el rango de fechas sin huecos.
"""
function validar_datos_entrada(datos_bundle::NamedTuple, tecnologias::Dict)
    # Validar demanda
    columnas_requeridas = [:anio, :estacion, :hora, :vector, :valor]
    columnas_presentes = Symbol.(names(datos_bundle.demanda))
    columnas_faltantes = setdiff(columnas_requeridas, columnas_presentes)
    
    if !isempty(columnas_faltantes)
        error("Faltan columnas en datos_bundle.demanda: $columnas_faltantes")
    end
    
    # Validar que no haya huecos en las fechas
    anios = sort(unique(datos_bundle.demanda.anio))
    for i in 2:length(anios)
        if anios[i] - anios[i-1] > 1
            @warn "Hueco detectado en años: $(anios[i-1]) a $(anios[i])"
        end
    end
    
    # Validar tecnologías
    if isempty(tecnologias)
        error("El diccionario de tecnologías está vacío")
    end
    
    # Validar que todas las tecnologías tengan campos requeridos
    for (nombre, tec) in tecnologias
        if !hasproperty(tec, :vector) || !hasproperty(tec, :costo_om_variable)
            error("Tecnología '$nombre' no tiene campos requeridos")
        end
    end
end

"""
    extraer_politicas_anio(politicas_df::DataFrame, anio::Int)

Extrae precio co2 y subvención máxima para un año específico.

# Retorna
NamedTuple con campos `precio_co2` y `subv_max_anual`
"""
function extraer_politicas_anio(politicas::DataFrame, anio::Int)
    filas = [ row for row in eachrow(politicas) if row.anio == anio ]
    if isempty(filas)
        error("No se encontró política para el año $anio")
    end
    
    return (
        precio_co2 = Float64(filas[1].precio_co2),
        subv_max_anual = Float64(filas[1].subv_max_anual)
    )
end

"""
    agregar_demanda_mensual(datos_bundle::NamedTuple, fechas::Vector{Date})

Agrega la demanda de todas las fechas especificadas y genera perfil horario repetido.

# Retorna
DataFrame con demanda horaria agregada para horizonte de 60 meses
"""
function agregar_demanda_mensual(datos_bundle::NamedTuple, fechas::Vector{Date})
    tmp_dia = DataFrame()
    for f in fechas
        append!(tmp_dia, Datos.demanda_dia(datos_bundle, f))
    end
    
    demanda_hora_mes = combine(groupby(tmp_dia, :hora),
                              :valor => sum => :demanda_mw)
    
    # Repetir perfil para horizonte de 60 meses
    return vcat([demanda_hora_mes for _ in 1:60]...)
end

function generar_perfil_aleatorizado(n_horas::Int)
    #μ, σ = 0.217, 0.5
    #dist = Normal(μ, σ)
    #perfil_dia = rand(dist, 24)
    #perfil_dia = clamp.(perfil_dia, 0.0, 1.0)
    #perfil_dia ./= maximum(perfil_dia)
    
    #exp_eolica = Dict{Int, Float64}()
    #for h in 0:(n_horas-1)
    #    exp_eolica[h] = perfil_dia[mod1(h, 24)]
    #end
    
    #return exp_eolica

    # Parámetros Beta para media ≈ 0.6 y varianza media - baja (aemet da el equivalente al 60% interpretándolo)
    α = 9.0
    β = 6.0  # Para que α/(α+β) ≈ 0.6
    
    dist_beta = Beta(α, β)
    
    # Generar perfil base de 24 horas
    perfil_dia = rand(dist_beta, 24)
    
    # Opcional: Añadir algo de autocorrelación temporal (suavizado simple)
    # Esto hace que los valores consecutivos sean más parecidos
    for i in 2:24
        perfil_dia[i] = 0.7 * perfil_dia[i] + 0.3 * perfil_dia[i-1]
        perfil_dia[i] = clamp(perfil_dia[i], 0.0, 1.0)  # Por seguridad
    end
        
    # Crear diccionario para n_horas
    exp_perfil = Dict{Int, Float64}()
    for h in 0:(n_horas-1)
        exp_perfil[h] = perfil_dia[mod1(h, 24)]
    end
        
    return exp_perfil
end



"""
    generar_perfil_eolico_aleatorio(n_horas::Int)

Genera un perfil de disponibilidad eólica aleatorio normalizado.

# Argumentos
- `n_horas::Int`: Número de horas para generar el perfil

# Retorna
Dict{Int,Float64} con factores de disponibilidad por hora [0,1]
"""
function generar_perfil_eolico_aleatorio(n_horas::Int)
    #μ, σ = 0.217, 0.5
    #dist = Normal(μ, σ)
    #perfil_dia = rand(dist, 24)
    #perfil_dia = clamp.(perfil_dia, 0.0, 1.0)
    #perfil_dia ./= maximum(perfil_dia)
    
    #exp_eolica = Dict{Int, Float64}()
    #for h in 0:(n_horas-1)
    #    exp_eolica[h] = perfil_dia[mod1(h, 24)]
    #end
    
    #return exp_eolica

    # Parámetros Beta para media ≈ 0.217 y varianza moderada
    α = 2.0
    β = 7.22  # Para que α/(α+β) ≈ 0.217
    
    dist_beta = Beta(α, β)
    
    # Generar perfil base de 24 horas
    perfil_dia = rand(dist_beta, 24)
    
    # Opcional: Añadir algo de autocorrelación temporal (suavizado simple)
    # Esto hace que los valores consecutivos sean más parecidos
    for i in 2:24
        perfil_dia[i] = 0.7 * perfil_dia[i] + 0.3 * perfil_dia[i-1]
        perfil_dia[i] = clamp(perfil_dia[i], 0.0, 1.0)  # Por seguridad
    end
        
    # Crear diccionario para n_horas
    exp_eolica = Dict{Int, Float64}()
    for h in 0:(n_horas-1)
        exp_eolica[h] = perfil_dia[mod1(h, 24)]
    end
        
    return exp_eolica
end

"""
    procesar_mercado_diario_integrado!(fecha, escenario, agentes, datos_bundle, tecnologias, precio_co2, produccion_anual)

Procesa el mercado diario para todos los vectores energéticos con recálculo iterativo.

# Retorna
Dict{Tuple{Int,String},Float64} con precios por (hora, vector)
"""
function procesar_mercado_diario_integrado!(
    fecha, escenario, agentes, datos_bundle, tecnologias, precio_co2_anio,
    produccion_anual, produccion_mensual, subv_unit_mes, mes, capacidad_disponible   
)
    # ———————————————————————————————
    # 0. Preparar estructuras de autoconsumo
    # autoconsumo_por_hora[vector][hora] = mwh consumidos
    autoconsumo_por_hora = Dict{String,Dict{Int,Float64}}(
        vector => Dict{Int,Float64}() for vector in obtener_vectores_ordenados(tecnologias)
    )
    produccion_dia = Dict{Int,Dict{String,Float64}}()
    # 1. Obtener demanda de todos los vectores
    # initial_daily_demands: La demanda base para el día, no se modifica dentro del bucle while.
    initial_daily_demands = Datos.demanda_dia_todos_vectores(datos_bundle, fecha)
    # current_dispatch_demands: Demanda usada para el despacho, se actualiza en cada iteración del while.
    # Para vectores no eléctricos, se basa en initial_daily_demands + una pasada de derivación.
    # Para electricidad, acumula nueva_demanda_extra.
    current_dispatch_demands = deepcopy(initial_daily_demands)

    # Obtener estación meteorológica
    estacion_fecha = obtener_estacion_fecha(datos_bundle, fecha)
    
    #precio_co2_anio_aux = precio_co2_anio[year(fecha)]
    

    generar_perfil_constante(v::Float64=1.0) = Dict(h => v for h in 0:23)
    # ── Construir perfiles aleatorios del día
    perf_al = Dict{String,Dict{Int,Float64}}()
    for tec in keys(tecnologias)
        if tec == "eolica"
            perf_al[tec] = generar_perfil_eolico_aleatorio(24)   # le asigno esta porque la otra quitaba mucha cantidad de potencia
        elseif tec in ("solar_pv", "solar_termica")
            #perf_al[tec] = generar_perfil_aleatorizado(24)
            perf_al[tec] = generar_perfil_constante()
        else
            perf_al[tec] = generar_perfil_constante()
        end
    end
    
    capacidad_disponible_hoy = Dict{Tuple{Int,Symbol,Int},Float64}()
    for (ag_id, ag) in agentes, t in keys(tecnologias), h in 0:23
        if !haskey(capacidad_disponible, (ag_id, Symbol(t), h))
            capacidad_disponible_hoy[(ag_id, Symbol(t), h)] = 0.0
        else
            capacidad_disponible_hoy[(ag_id, Symbol(t), h)] = capacidad_disponible[(ag_id, Symbol(t), h)] * perf_al[t][h]
        end
    end


    #for (ag_id, ag) in agentes, t in keys(tecnologias), h in 0:23
    #    if !haskey(capacidad_disponible, (ag_id, Symbol(t), h))
    #        capacidad_disponible_hoy[(ag_id, Symbol(t), h)] = 0.0
    #    else
    #        capacidad_disponible_hoy[(ag_id, Symbol(t), h)] *= perf_al[t][h]
    #    end
    #end

    # ------------------------------------------------------------------
    # 1. Curva de aprendizaje sobre el coste marginal variable
    # ------------------------------------------------------------------
    for tec in values(tecnologias)
        Tecnologias.actualizar_costo_marginal_aprendizaje!(tec)
    end

    # Coste variable €/MWh tomado ya con la reducción por aprendizaje
    costo_marginal = Dict{Tuple{Int,Symbol},Float64}()

    for (ag_id, ag) in agentes
        for (tec, _mw) in ag.capacidades               # basta con que el agente tenga la tecnología
            if haskey(tecnologias, tec)
                # ------------------------------------------
                # Coste variable total = O&M + CO₂
                # ------------------------------------------
                cv_base = tecnologias[tec].costo_om_variable      # €/MWh
                emis    = tecnologias[tec].emisiones_unitarias     # t CO₂ / MWh
                costo_co2 = emis * precio_co2_anio                # €/MWh
                cv_total = cv_base + costo_co2
            else
                @warn "Tecnología $tec no encontrada en diccionario tecnologías; cv=0"
                cv_total = 0.0
            end
            costo_marginal[(ag_id, Symbol(tec))] = cv_total
        end
    end

    for (ag_id, _) in agentes        # recorre todos los agentes
        for t in keys(tecnologias)   # recorre todas las tecnologías
            clave = (ag_id, Symbol(t))
            if !haskey(costo_marginal, clave)
                # toma el coste OM de la tecnología como valor por defecto
                costo_marginal[clave] = tecnologias[t].costo_om_variable
            end
        end
    end


    demanda_dia = DataFrame(hora = Int[], vector = String[], valor = Float64[])
    for h in 0:23
        for (vec, mwh) in current_dispatch_demands[h]
            push!(demanda_dia, (h, vec, mwh))
        end
    end



    # Antes de aplanar, rellena con zeros
    all_vectors = unique([tec.vector for tec in values(tecnologias)])
    for h in keys(current_dispatch_demands)
        for vec in all_vectors
            current_dispatch_demands[h][vec] = 
                get(current_dispatch_demands[h], vec, 0.0)
        end
    end

    # ─────────────────────────────────────────────────────────────
    # Demanda plana {(:vector, hora) => MWh} que exige el solver
    # ─────────────────────────────────────────────────────────────
    demand_tuple = Dict{Tuple{Symbol,Int},Float64}()
    for (h, dict_vec) in current_dispatch_demands
        for (vec, mwh) in dict_vec
            demand_tuple[(Symbol(vec), h)] = mwh
        end
    end

    # ────────────── WARM-START ──────────────
    if month(fecha) == 1 && day(fecha) == 1
        x0_prev = nothing # DESTRUYO LA SEMILLA PARA INICIO DE AÑO. COMO CAMBIA TANTO LA DEMANDA Y, ADEMÁS, EN EL PRIMER AÑO
    else                    # ENTRA EN JUEGO LA PRODUCCIÓN QUE TENÍA RETRASO DE 12 MESES, LA SEMILLA LÍA AL MODELO, ES MUY DIFF.
        x0_prev = get(LAST_LCP_SOLUTION, fecha - Day(1), nothing)
    end

    println("📅 Día: ", day(fecha), ", Mes: ", month(fecha), ", Año: ", year(fecha)) #REVISAR, INSERTAR SIMBOLITO, ESTÁ MUY TRISTE
    resultado = EquilibrioMultivector.subastar_multivector!(
        fecha,
        agentes, tecnologias,
        demand_tuple,
        capacidad_disponible_hoy,
        costo_marginal;
        x0 = x0_prev)

    # Guardamos la solución para la próxima jornada
    LAST_LCP_SOLUTION[fecha] = resultado.x0_next

    precios_todos_vectores = resultado.precios_marginales
    precios_convertidos = Dict{Tuple{Int,String},Float64}()
    for ((vec_sym, h), precio) in precios_todos_vectores
        precios_convertidos[(h, String(vec_sym))] = precio
    end # PARA RESPETAR COMO ESTABA DENOMINADO EN EL OTRO ESPACIO

    
    # ───────────────────────────────────────────────────────────────────
    # Aplanar despacho (agente, tecnología, vector, hora)
    # a producción por agente-tecnología-hora
    dispatch_flat = resultado.despacho                  # Dict{Tuple{Int,String,String,Int},Float64}
    produccion_dia = Dict{Int,Dict{String,Float64}}()

    # ── 4-bis.  Registrar finanzas de cada operación diaria ───────────────
    subv_del_mes = get(subv_unit_mes, month(fecha), Dict{String,Float64}())
    for ((ag_id, tec, vector, h), mwh) in dispatch_flat
        mwh < 1e-6 && continue
        ag  = agentes[ag_id]
        t   = tecnologias[tec]
        precio_venta = precios_convertidos[(h, vector)]
        costo_om_var = get(ag.var_costs, tec, t.costo_om_variable)
        costo_emis   = t.emisiones_unitarias * precio_co2_anio
        subv_unit    = get(subv_del_mes, tec, 0.0)
        Agentes.registrar_operacion_diaria!(
            ag,
            mwh,
            precio_venta,
            costo_om_var,
            costo_emis,
            subv_unit,
        )
    end
    # ── 4-ter. Acumular producción diaria por agente y tecnología ────
    for ((ag_id, tec, _vec, h), mwh) in dispatch_flat
        mwh < 1e-6 && continue
        # Inicializa el sub-dict si hace falta
        if !haskey(produccion_dia, ag_id)
            produccion_dia[ag_id] = Dict{String,Float64}()
        end
        # Suma la producción en MWh
        produccion_dia[ag_id][tec] =
            get(produccion_dia[ag_id], tec, 0.0) + mwh
    end
    # ───────────────────────────────────────────────────────────────────
    deficits_totales       = resultado.deficits
    if day(fecha) % 32 == 1
        println("Deficits totales: ", deficits_totales)
        #println("Precios totales: ", precios_todos_vectores)
        for (ag_id, ag) in agentes, tec_name in keys(tecnologias), h in 0:23
            clave = (ag_id, tec_name, h)
            if haskey(produccion_dia, clave) && produccion_dia[clave] > 0
                println("  Agente $(ag_id), tecnología $(tec_name), hora $(h): ",
                        produccion_dia[clave], " MWh")
            end
        end
    end

    # Actualizar producción diaria de agentes con el estado final de produccion_dia
    actualizar_produccion_diaria_agentes!(agentes, produccion_dia, fecha)
    

    # NUEVO: acumular en los diccionarios globales para informes
    acumular_produccion!(produccion_anual, produccion_mensual,
                         fecha, produccion_dia)
    # Guardar autoconsumo por fecha (estado final de autoconsumo_por_hora)
    #FULL_AUTOCONSUMO[fecha] = copy(autoconsumo_por_hora)
    #estado_iteracion[:produccion_acumulada] = produccion_dia
    #estado_iteracion[:precios_convergidos] = !pendiente_recalculo

    return (
        precios  = precios_convertidos, # Precios finales de esta iteración
        autoconsumo = copy(autoconsumo_por_hora), # Autoconsumo final
        deficits = deficits_totales
    )
end

"""
    obtener_estacion_fecha(datos_bundle::NamedTuple, fecha::Date)

Obtiene la estación meteorológica para una fecha dada.

# Retorna
String con el nombre de la estación ("invierno", "primavera", "verano", "otono") o "promedio"
"""
function obtener_estacion_fecha(datos_bundle::NamedTuple, fecha::Date)
    # Primero intentar con Dict (caso más común)
    if isa(datos_bundle.cal_estaciones, Dict)
        return get(datos_bundle.cal_estaciones, fecha, "promedio")
    end
    
    # Si es DataFrame, buscar la fila correspondiente
    if isa(datos_bundle.cal_estaciones, DataFrame)
        fila = filter(row -> row.fecha == fecha, datos_bundle.cal_estaciones)
        return isempty(fila) ? "promedio" : String(fila[1, :estacion])
    end
    
    # Fallback
    return "promedio"
end

"""
    calcular_precios_todos_vectores(precios_dict::Dict, tecnologias::Dict, agentes::Dict)

Calcula los precios de todos los vectores energéticos basado en costos marginales.

# Retorna
Dict{Tuple{Int,String},Float64} con precios por (hora, vector)
"""
function calcular_precios_todos_vectores(precios_dict::Dict, tecnologias::Dict, agentes::Dict)
    precios_todos = Dict{Tuple{Int,String},Float64}()
    
    for hora in 0:23
        # 1️⃣ Precio base (electricidad) SIN cambiar la variable del diccionario
        precio_elec_hora = precios_dict[(hora, "electricidad")]
        precios_todos[(hora, "electricidad")] = precio_elec_hora

        # 2️⃣ Resto de vectores
        for vector in obtener_vectores_ordenados(tecnologias)
            if vector != "electricidad"
                precios_todos[(hora, vector)] =
                    calcular_precio_vector_con_dependencias(
                        vector, hora, precio_elec_hora,
                        tecnologias, precios_todos, agentes)
            end
        end
    end
    
    return precios_todos
end

"""
    despachar_vector_no_electrico!(vector, demanda_mwh, hora, tecnologias, agentes, precios, produccion_dia, produccion_anual, fecha, capacidad_comprometida)

Despacha la demanda de un vector no eléctrico usando orden de mérito económico.

# Argumentos
- `vector::String`: Vector energético a despachar
- `demanda_mwh::Float64`: Demanda en mwh
- `hora::Int`: Hora del día (0-23)
- `tecnologias::Dict`: Diccionario de tecnologías
- `agentes::Dict`: Diccionario de agentes
- `precios::Dict`: Precios por (hora, vector)
- `produccion_dia::Dict`: Acumulador de producción diaria
- `produccion_anual::Dict`: Acumulador de producción anual
- `fecha::Date`: Fecha del despacho

# Retorna
Float64 con el consumo eléctrico adicional generado (mwh)
"""
function despachar_vector_no_electrico!(
    vector, demanda_mwh, hora, tecnologias, agentes, precios,
    produccion_dia, produccion_anual, fecha,
    autoconsumo_por_hora::Dict{String,Dict{Int,Float64}},
    precio_co2::Float64,
    capacidad_comprometida::Dict{Tuple{Int,String,Int},Float64}
)::Tuple{Dict{String,Float64}, Float64}
    # Encontrar productores del vector
    productores = construir_lista_productores(vector, hora, tecnologias, agentes, precios, precio_co2)
    
    # Ordenar por costo
    sort!(productores, by = x -> x[1])
    
    # Despachar en orden de mérito
    demanda_restante = demanda_mwh
    consumo_elec_total = 0.0
    consumos_por_vector = Dict{String,Float64}()

    for (costo, ag_id, tec_name, cap_disp, tec) in productores
        if demanda_restante < 1e-6
            break
        end

        clave        = (ag_id, tec_name, hora)
        cap_ya_usada = get(capacidad_comprometida, clave, 0.0)
        cap_disponible_real = max(0.0, cap_disp - cap_ya_usada)
        despachado = min(demanda_restante, cap_disponible_real)
        capacidad_comprometida[clave] = cap_ya_usada + despachado
        if !haskey(produccion_dia, ag_id)
            produccion_dia[ag_id] = Dict{String,Float64}()
        end
        produccion_dia[ag_id][tec_name] =
            get(produccion_dia[ag_id], tec_name, 0.0) + despachado

            anio = year(fecha)
            if !haskey(produccion_anual, anio)
                produccion_anual[anio] = Dict{String,Float64}()
            end
            produccion_anual[anio][tec_name] =
                get(produccion_anual[anio], tec_name, 0.0) + despachado        

        if !isempty(tec.vector_consumido) && tec.ratio_energia > 0
            consumo_vector = despachado / tec.ratio_energia
            consumos_por_vector[tec.vector_consumido] = 
                get(consumos_por_vector, tec.vector_consumido, 0.0) + consumo_vector
        end

        demanda_restante -= despachado
        
        # Asignar producción a agentes proporcionalmente
        #consumo_elec_tec = asignar_produccion_agentes!(
        #    despachado, tec_name, cap_disp, tec, vector,
        #    agentes, precios, hora, produccion_dia, fecha,
        #    autoconsumo_por_hora
        #)
        
        #consumo_elec_total += consumo_elec_tec
        # producción acumulada del agente y tecnología
        consumo_elec_tec = asignar_produccion_agentes!(
            despachado, tec_name, cap_disp, tec, vector,
            agentes, precios, hora, produccion_dia, fecha,
            autoconsumo_por_hora
        )
        consumo_elec_total += consumo_elec_tec
        # — acumuladores anuales y mensuales —
        anio = year(fecha)
        produccion_anual[anio][tec_name] =
            get(produccion_anual[anio], tec_name, 0.0) + despachado
        pm = get!(produccion_mensual[anio], month(fecha),
                  Dict{String,Float64}())
        pm[tec_name] = get(pm, tec_name, 0.0) + despachado

    end
    
    deficit = max(0.0, demanda_restante)


    if deficit > 1e-3
        # REVISAR LO QUITO AHORA PARA AGILIZAR println("    ⚠️  Demanda no cubierta para $vector H$hora: $(round(deficit, digits=2)) mwh")
    end
    


    return consumos_por_vector, deficit
end

"""
    construir_lista_productores(vector, hora, tecnologias, agentes, precios)

Construye lista de productores disponibles para un vector con sus costos marginales.
"""
function construir_lista_productores(
    vector::String, hora::Int,
    tecnologias::Dict, agentes::Dict, precios::Dict,
    precio_co2::Float64
)
    productores = []
    
    for (tec_name, tec) in tecnologias
        if tec.vector == vector
            # Calcular costo marginal
            costo_base = tec.costo_om_variable #€/mwh
            costo_input = 0.0
            if !isempty(tec.vector_consumido) && !ismissing(tec.ratio_energia) &&
               tec.ratio_energia > 1e-6
                precio_input = get(precios, (hora, tec.vector_consumido), 0.0)
                costo_input = precio_input / tec.ratio_energia
            end

            # emisiones_unitarias está en toneladasco2/mwh; precio_co2 en €/tco2
            emis_t = tec.emisiones_unitarias # → tco2/mwh
            costo_co2   = precio_co2 * emis_t
            costo_total = costo_base + costo_input + costo_co2
            
            
            # Capacidad disponible
            cap_disponible = sum(get(ag.capacidades, tec_name, 0.0) for ag in values(agentes))
            
            if cap_disponible > 1e-6
                for ag in values(agentes)
                    cap_ag = get(ag.capacidades, tec_name, 0.0)
                    cap_ag > 1e-6 &&
                        push!(productores, (costo_total, ag.id, tec_name, cap_ag, tec))
                end
            end
            if vector == "gas" && cap_disponible > 1e-6
                #println("      🔍 Productor de gas encontrado: $tec_name con capacidad $cap_disponible MW")
            end
        end
    end

    return productores
end

"""
    asignar_produccion_agentes!(despachado, tec_name, cap_disp, tec, vector, agentes, precios, hora, produccion_dia, fecha)

Asigna producción despachada a los agentes proporcionalmente y actualiza sus finanzas.

# Retorna
Float64 con el consumo eléctrico si la tecnología consume electricidad
"""
function asignar_produccion_agentes!(
    despachado, tec_name, cap_disp, tec, vector, agentes, precios, hora,
    produccion_dia, fecha,
    autoconsumo_por_hora::Dict{String,Dict{Int,Float64}}
)::Float64
    consumo_elec_total = 0.0
    
    for ag in values(agentes)
        cap_ag = get(ag.capacidades, tec_name, 0.0)
        if cap_ag > 0
            proporcion = cap_ag / cap_disp
            prod_ag = despachado * proporcion
            
            # Actualizar finanzas del agente
            precio_venta = get(precios, (hora, vector), 0.0)
            ag.revenues += prod_ag * precio_venta
            ag.cost += prod_ag * tec.costo_om_variable
            # Coste energético del insumo
            if !isempty(tec.vector_consumido) && tec.ratio_energia > 0
                precio_in = get(precios, (hora, tec.vector_consumido), 0.0)    # €/mwh-input
                coste_eur = (prod_ag * precio_in / tec.ratio_energia)
                ag.cost += coste_eur
            end



            # Actualizar producción anual del agente (en mwh)
            if !hasproperty(ag, :produccion_anual)
                ag.produccion_anual = Dict{String,Float64}()
            end
            ag.produccion_anual[tec_name] = get(ag.produccion_anual, tec_name, 0.0) + prod_ag
            
            # Registrar producción diaria
            if !haskey(produccion_dia, ag.id)
                produccion_dia[ag.id] = Dict{String,Float64}()
            end
            produccion_dia[ag.id][tec_name] = get(produccion_dia[ag.id], tec_name, 0.0) + prod_ag
        end
    end
    
    # ————— Registrar autoconsumo del vector_consumido (eléctrico u otro)
    consumo_vector = 0.0
    if !isempty(tec.vector_consumido) && tec.ratio_energia > 0
        # ratio_energia ahora se interpreta como ratio de insumo genérico
        consumo_vector = despachado / tec.ratio_energia
        # actualizar autoconsumo_por_hora
        apre = autoconsumo_por_hora[tec.vector_consumido]
        apre[hora] = get(apre, hora, 0.0) + consumo_vector
    end
    
    # Devolver consumo para su registro en la fase de recálculo
    return consumo_vector
end

"""
    actualizar_produccion_diaria_agentes!(agentes, produccion_dia, fecha)

Actualiza la producción diaria registrada en cada agente.
"""
function actualizar_produccion_diaria_agentes!(
        agentes::Dict, produccion_dia::Dict, fecha::Date)

    anio = year(fecha)
    mes  = month(fecha)
    for (ag_id, prod_tec) in produccion_dia
        haskey(agentes, ag_id) || continue
        agente = agentes[ag_id]

        # ── 1. Guardar desglose diario ───────────────────────────────────
        if !hasproperty(agente, :produccion_diaria)
            agente.produccion_diaria = Dict{Date,Dict{String,Float64}}()
        end
        agente.produccion_diaria[fecha] = prod_tec

        # ── 2. Acumular en anual y mensual ───────────────────────────────
        if !hasproperty(agente, :produccion_anual)
            agente.produccion_anual = Dict{String,Float64}()
        end
        if !hasproperty(agente, :produccion_mensual)
            agente.produccion_mensual =
                Dict{Int,Dict{Int,Dict{String,Float64}}}()
        end
        # Diccionario mensual ↴
        pm_ag   = get!(agente.produccion_mensual, anio,
                        Dict{Int,Dict{String,Float64}}())
        pm_mes  = get!(pm_ag, mes, Dict{String,Float64}())

        # Recorremos las tecnologías producidas en el día
        for (tec, mwh) in prod_tec
            # Acumulador anual
            agente.produccion_anual[tec] =
                get(agente.produccion_anual, tec, 0.0) + mwh
            # Acumulador mensual
            pm_mes[tec] = get(pm_mes, tec, 0.0) + mwh
        end
    end
end


"""
    obtener_vectores_ordenados(tecnologias::Dict{String,Tecnologia})

Construye un orden para competir los vectores:
1. Agrupa en SCCs (ciclos) → cada componente es un bloque.
2. Topological sort de esos bloques.
3. Dentro de cada bloque, orden alfabético de vectores.
"""
function obtener_vectores_ordenados(tecnologias::Dict{String,Tecnologia})
    # 1) Recopilar nodos y aristas
    vects = Set{String}()
    aristas = Vector{Tuple{String,String}}()
    for tec in values(tecnologias)
        vec = tec.vector
        isempty(vec) && continue
        push!(vects, vec)
        vcons = tec.vector_consumido
        isempty(vcons) && continue
        push!(vects, vcons)
        push!(aristas, (vcons, vec))
    end

    # 2) Índices para Graphs
    lista = collect(vects)
    idx = Dict(v=>i for (i,v) in enumerate(lista))
    g = DiGraph(length(lista))
    for (u,v) in aristas
        add_edge!(g, idx[u], idx[v])
    end

    # 3) SCCs y mapeo nodo→comp
    comps = strongly_connected_components(g)
    comp_of = Dict{Int,Int}()
    for (cid, comp) in enumerate(comps)
        for v in comp
            comp_of[v] = cid
        end
    end

    # 4) Construir grafo de componentes
    ncomp = length(comps)
    gc = DiGraph(ncomp)
    for (u,v) in aristas
        cu, cv = comp_of[idx[u]], comp_of[idx[v]]
        cu != cv && add_edge!(gc, cu, cv)
    end

    # 5) Topological sort de componentes
    orden_comp = topological_sort(gc)

    # 6) Desplegar vectores: para cada comp en orden, sus vectores α-ordenados
    resultado = String[]
    for cid in orden_comp
        # Extrae los nombres de los nodos de la componente cid
        nombres = [ lista[v] for v in comps[cid] ]
        # Orden alfabético dentro de la componente
        append!(resultado, sort(nombres))
    end

    return resultado
end



"""
    calcular_precio_vector_con_dependencias(vector, hora, precio_input, tecnologias, precios_conocidos, agentes)

Calcula el precio de un vector considerando sus dependencias de otros vectores.

# Retorna
Float64 con el precio calculado en €/mwh
"""
function calcular_precio_vector_con_dependencias(
    vector::String, hora::Int, precio_input::Float64,
    tecnologias::Dict, precios_conocidos::Dict, agentes::Dict
)
    # Encontrar tecnología más barata que produce este vector
    costo_minimo = Inf
    
    for (tec_name, tec) in tecnologias
        if tec.vector == vector
            # Verificar si algún agente tiene capacidad
            cap_total = sum(get(ag.capacidades, tec_name, 0.0) for ag in values(agentes))
            if cap_total < 1
                continue
            end
            
            costo = tec.costo_om_variable
            
            if !isempty(tec.vector_consumido)
                if tec.vector_consumido == "electricidad"
                    precio_input = precio_input
                else
                    precio_input = get(precios_conocidos, (hora, tec.vector_consumido), precio_input * 0.8)
                end
                
                if tec.ratio_energia > 0.01
                    costo += precio_input / tec.ratio_energia
                end
            end
            
            costo_minimo = min(costo_minimo, costo)
        end
    end
    
    # Si no hay productores, usar precio base
    if costo_minimo == Inf
        return precio_input * 0.8
    end
    
    # Añadir markup del 10%
    return costo_minimo * 1.1
end

"""
    construir_snapshot_finanzas(agentes, anio)

Construye un snapshot de las finanzas de todos los agentes para un año.
"""
function construir_snapshot_finanzas(agentes::Dict, anio::Int)
    return [
        (
            agent_id = ag.id,
            anio = anio,
            empresa = ag.name,
            cash_eur = ag.cash,
            debt_eur = ag.debt,
            revenues_eur = ag.revenues,
            cost_eur = ag.cost,
            profit_eur = ag.profit,
            ingresos = ag.revenues,
            costo_operacion = ag.cost,
            beneficio_neto = ag.profit
        )
        for ag in values(agentes)
    ]
end

"""
    agregar_capacidades_sistema(agentes, tecnologias)

Agrega las capacidades totales del sistema por tecnología.
"""
function agregar_capacidades_sistema(agentes::Dict, tecnologias::Dict)
    return Dict(
        tecnologia => sum(ag.capacidades[tecnologia] for ag in values(agentes) 
                        if haskey(ag.capacidades, tecnologia); init=0.0)
        for tecnologia in keys(tecnologias)
    )
end

"""
    construir_resultados_simulacion(escenario, agentes, capacidades, produccion, finanzas, precios_hist, datos_bundle)

Construye la estructura final de resultados de la simulación.
"""
function construir_resultados_simulacion(
    escenario, agentes, capacidades_finales,
    produccion_anual, finanzas_anuales, precios_historicos,
    produccion_mensual::Dict, deficits_mensual::Dict, datos_bundle::NamedTuple
)
    # Construir DataFrames de resultados
    df_capacidad = construir_df_capacidades(capacidades_finales)
    df_produccion = construir_df_produccion(produccion_anual)
    
    # --- construir DataFrame de producción mensual ---
    df_produccion_mensual = DataFrame(
      anio       = Int[],
      mes        = Int[],
      tecnologia = String[],
      valor      = Float64[]
    )
    for (anio, meses) in produccion_mensual
      for (mes, dict_tec) in meses
        for (tec, mwh) in dict_tec
          push!(df_produccion_mensual, (anio, mes, tec, mwh))
        end
      end
    end

    # --- construir DataFrame de déficit mensual ---
    df_deficit_mensual = DataFrame(
      anio   = Int[],
      mes    = Int[],
      vector = String[],
      valor  = Float64[]
    )
    for (anio, meses) in deficits_mensual
      for (mes, dict_vec) in meses
        for (vec, mwh) in dict_vec
          push!(df_deficit_mensual, (anio, mes, vec, mwh))
        end
      end
    end

    # Construir diccionario de precios para compatibilidad
    precios_sistema = Dict{Tuple{Int,String,String,Int}, Float64}()
    for (fecha, precios_dia) in precios_historicos
        anio_precio = year(fecha)
        estacion_precio = obtener_estacion_fecha(datos_bundle, fecha)
        for ((hora, vector), precio) in precios_dia
            precios_sistema[(anio_precio, vector, estacion_precio, hora)] = precio
        end
    end
    
    beneficio_total_sistema = sum(ag.profit for ag in values(agentes)) / 1e6
    
    return (
        escenario = escenario.nombre,
        opt_status = MOI.OPTIMAL,
        opt_objective = beneficio_total_sistema, # MILLONES €
        opt_resultados = Dict(
            :capacidad  => df_capacidad,
            :produccion => df_produccion,
            :produccion_mensual => df_produccion_mensual,
            :deficit_mensual => df_deficit_mensual,
        ),
        agentes = agentes,
        finanzas_anuales = finanzas_anuales,
        precios = precios_sistema,
        error_info = nothing
    )
end

# Helper functions
function get_precio_co2_fecha(escenario, fecha::Date)::Float64
    anio = year(fecha)
    mes = month(fecha)
    return get(escenario.precio_co2, (anio, mes), 0.0)
end

function construir_df_capacidades(caps::Dict{Int,Dict{String,Float64}})
    filas = [(anio, tecnologia, valor, 0.0) for (anio, dict_tec) in caps for (tecnologia, valor) in dict_tec if valor > 1e-6]
    return DataFrame(
        anio = [f[1] for f in filas],
        tecnologia = [f[2] for f in filas],
        cap_total = [f[3] for f in filas],
        cap_nueva = [f[4] for f in filas]
    )
end

function construir_df_produccion(prod::Dict{Int,Dict{String,Float64}})
    filas = [(anio, tecnologia, "anual", 0, valor) for (anio, dict_tec) in prod for (tecnologia, valor) in dict_tec if valor > 1e-6]
    return DataFrame(
        anio = [f[1] for f in filas],
        tecnologia = [f[2] for f in filas],
        estacion = [f[3] for f in filas],
        hora = [f[4] for f in filas],
        valor = [f[5] for f in filas]
    )
end



# Usar la lógica de conversión que ya existe en OptimizacionSistema.jl
function calcular_precio_usando_conversion_existente(vector::String, precio_base::Float64, agente::Agentes.EnergyAgent, tecnologias::Dict)
    # Buscar tecnologías que producen este vector en el agente
    precio_minimo = precio_base
    
    for tec_name in keys(agente.capacidades)
        if haskey(tecnologias, tec_name)
            tec = tecnologias[tec_name]
            if tec.vector == vector
                # Usar el costo variable del agente (ya calculado)
                costo_var = get(agente.var_costs, tec_name, 0.0)
                precio_ajustado = costo_var  #€/mwh
                
                # Si consume otro vector, aplicar factor de conversión
                if !isempty(tec.vector_consumido)
                    precio_ajustado = precio_ajustado / max(tec.ratio_energia, 0.1)
                end
                
                precio_minimo = min(precio_minimo, precio_ajustado)
            end
        end
    end
    
    return max(precio_minimo, precio_base * 0.3)
end

function calcular_precio_por_vector_secuencial(target_vector::String, precio_electricidad::Float64, tecnologias::Dict, all_current_vector_prices::Dict{String, Float64})
    # Encontrar tecnologías que producen este vector
    productoras = [ pair for pair in collect(tecnologias) if pair.second.vector == target_vector ]
    
    if isempty(productoras)
        # Fallback if no direct producers found for the target_vector
        # This might need adjustment based on broader system design,
        # but for now, mirrors the old behavior relative to electricity price.
        return precio_electricidad * 0.8 
    end
    
    costos_marginales = Dict{String, Float64}()
    for (tec_nombre, tec) in productoras # Changed to use tec_nombre for logging
        costo_produccion_propio = tec.costo_om_variable # €/mwh of output vector
        costo_input_energia = 0.0

        if !isempty(tec.vector_consumido)
            input_vector_price = 0.0
            if tec.vector_consumido == "electricidad"
                input_vector_price = precio_electricidad
            elseif haskey(all_current_vector_prices, tec.vector_consumido)
                input_vector_price = all_current_vector_prices[tec.vector_consumido]
            else
                # Consumed vector's price is not yet known for this hour.
                println("WARN: Price for consumed vector $(tec.vector_consumido) for technology $(tec_nombre) not found in current iteration. Skipping this technology for costing $(target_vector).")
                continue # Skips this 'tec' and goes to the next in the loop
            end

            if tec.ratio_energia > 1e-6 # Avoid division by zero or very small ratio
                costo_input_energia = input_vector_price / tec.ratio_energia
            else
                println("WARN: ratio_energia for technology $(tec_nombre) is invalid (<= 1e-6). Skipping this technology for costing $(target_vector).")
                continue # Skips this 'tec'
            end
        end
        costos_marginales[tec_nombre] = costo_produccion_propio + costo_input_energia
    end
    
    if isempty(costos_marginales)
        # This case occurs if all producers were skipped (e.g., due to missing input prices or invalid ratios)
        # Returning a high price or a specific error indicator might be better.
        # For now, retain existing fallback logic relative to electricity, though this might need review.
        # Consider what an appropriate fallback is if no producer can be costed.
        println("WARN: No valid producers found for vector $(target_vector) after costing. Returning fallback price.")
        return precio_electricidad * 0.8 # Fallback, same as when no producers
    end
    
    return minimum(costos_marginales)
end

"""
Calcula el precio de escasez para un vector con déficit.
Aplica una prima de escasez proporcional al déficit.
"""
function calcular_precio_escasez(precio_base::Float64, deficit::Float64, demanda_total::Float64)::Float64
    if demanda_total < 1e-6 || deficit < 1e-6
        return precio_base
    end
    

    # Calcular ratio de déficit
    ratio_deficit = 1 + 100 * deficit / (demanda_total - deficit + 0.01) #deficit / demanda_total
    if deficit != demanda_total

        # Prima de escasez: crece exponencialmente con el déficit
        # Con 10% déficit → 50% más caro
        # Con 50% déficit → 300% más caro
        # Con 100% déficit → 1000% más caro (10x)
        prima_escasez = precio_base * ratio_deficit # (exp(3 * ratio_deficit)) QUITO ESTO PORQUE VOY A PROBAR EL MULTIPLICADOR DEL VICTORIA 3
        
        # Precio mínimo de escasez: 200 €/MWh
        precio_escasez = min(max(precio_base + prima_escasez, 500.0),20000.0)
    
    else
        precio_escasez = 20000.0
    end


    # Cap máximo: 3000 €/MWh (evitar valores absurdos)
    return min(precio_escasez, 20000.0)
end


"""
Aplica penalizaciones económicas a los agentes por déficit de suministro.
Las penalizaciones se distribuyen proporcionalmente a la cuota de mercado teórica.
"""
function aplicar_penalizaciones_deficit!(
    agentes::Dict{Int,Agentes.EnergyAgent},
    deficits::Dict{Tuple{Int,String},Float64},
    fecha::Date,
    tecnologias::Dict{String,Tecnologia}
)
    if isempty(deficits)
        return
    end
    
    # Calcular penalizaciones por vector
    penalizaciones_por_agente = Dict{Int,Float64}()
    
    for ((hora, vector), deficit) in deficits
        # CAMBIO: Penalizar a TODOS los agentes proporcionalmente
        # a su participación en el mercado total
        
        capacidad_total_sistema = 0.0
        capacidades_por_agente = Dict{Int,Float64}()
        
        # Calcular capacidad total de TODOS los agentes (no solo del vector)
        for (ag_id, ag) in agentes
            cap_ag_total = sum(values(ag.capacidades); init=0.0)
            capacidades_por_agente[ag_id] = cap_ag_total
            capacidad_total_sistema += cap_ag_total
        end
        
        if capacidad_total_sistema < 1e-6
            # Sistema sin capacidad: distribuir equitativamente
            n_agentes = length(agentes)
            penalizacion_por_agente = (deficit )  / n_agentes # *1000 (Que al menos salga penalizado 1000€ el mw)
            
            for ag_id in keys(agentes)
                penalizaciones_por_agente[ag_id] = get(penalizaciones_por_agente, ag_id, 0.0) + penalizacion_por_agente
            end
        else
            # Distribuir proporcionalmente a la capacidad total del agente
            # (incentiva a todos a cubrir déficits)
            for (ag_id, cap_ag) in capacidades_por_agente
                proporcion = cap_ag / capacidad_total_sistema
                penalizacion = deficit * proporcion #*1000 #(Que al menos salga penalizado 1000€ el mw)
                penalizaciones_por_agente[ag_id] = get(penalizaciones_por_agente, ag_id, 0.0) + penalizacion
            end
        end
    end
    
    # Aplicar penalizaciones...
    for (ag_id, penalizacion) in penalizaciones_por_agente
        if penalizacion > 1e-6
            agentes[ag_id].cost += penalizacion
            #agentes[ag_id].profit -= penalizacion SE PENALIZABA 2 VECES
            agentes[ag_id].cash -= penalizacion * 0.05  # 5% se paga inmediatamente
            println("      ⚠️  Agente $(agentes[ag_id].name): penalización $(round((penalizacion/1e6), digits=2)) M€ por déficit")


            if !hasproperty(agentes[ag_id], :historial_penalizaciones) ||
                agentes[ag_id].historial_penalizaciones === nothing
                 agentes[ag_id].historial_penalizaciones =
                     Dict{Tuple{Int,Int},Float64}()
             end
            anio_pen = year(fecha); mes_pen = month(fecha)
            agentes[ag_id].historial_penalizaciones[(anio_pen, mes_pen)] =
                    get(agentes[ag_id].historial_penalizaciones,
                        (anio_pen, mes_pen), 0.0) + penalizacion
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
#  NUEVAS FUNCIONES  ➜  MILP rolling-horizon + fallback heurístico
# ──────────────────────────────────────────────────────────────────────────────
            
"""
invertir_mensualmente_con_milp!( … )
            
Ejecuta el MILP multi-agente con un horizonte rodante de 5 años y aplica  
únicamente la inversión del mes actual. Si HiGHS devuelve :Infeasible,  
timeout o cualquier error, cae al fallback `PlanificacionMes`.
"""
function invertir_mensualmente_con_milp!(
    agentes::Dict{Int,Agentes.EnergyAgent},
    tecs::Dict{String,Tecnologias.Tecnologia},
    datos_bundle::NamedTuple,
    esc::Escenario,
    anio::Int,
    mes::Int,
    precios_historicos::Dict{Date,Dict{Tuple{Int,String},Float64}},
    deficits_por_vector::Dict{String,Float64},
    precio_co2_anio::Float64,
    precio_medio_vector::Dict{String,Float64}
)
    # 1️⃣  Horizonte rolling de 5 años
    HORIZONTE = 5
    anio_fin  = min(anio + HORIZONTE - 1, 2050)

    datos_horizonte = (
        demanda             = filter(r -> anio <= r.anio <= anio_fin, datos_bundle.demanda),
        perfiles_generacion = datos_bundle.perfiles_generacion,
        politicas           = filter(r -> anio <= r.anio <= anio_fin, datos_bundle.politicas),
        tasa_descuento      = datos_bundle.tasa_descuento,
        cal_estaciones      = datos_bundle.cal_estaciones,
        potencial_renovable = datos_bundle.potencial_renovable,
        delay_construccion  = datos_bundle.delay_construccion,
        perfiles_gen        = datos_bundle.perfiles_gen,
        exp_eolica          = datos_bundle.exp_eolica,
        precio_co2_anio     = precio_co2_anio,
        precio_medio_vector = precio_medio_vector,
    )

    # ——— Precio esperado derivado de los históricos ya observados ———
    precios_mercado_horizonte = Dict{Tuple{Int,String,String,Int},Float64}()

    if !isempty(precios_historicos)
        # Media de los ÚLTIMOS 30 días para cada (vector, estación, hora)
        acumulados   = Dict{Tuple{String,String,Int},Tuple{Float64,Int}}()
        fechas_rec   = sort(collect(keys(precios_historicos)), rev = true)[1:min(30, length(precios_historicos))]

        for f in fechas_rec
            est = get(datos_bundle.cal_estaciones, f, "desconocido")
            for ((hora, vec), precio) in precios_historicos[f]
                k            = (vec, est, hora)
                suma, n      = get(acumulados, k, (0.0, 0))
                acumulados[k] = (suma + precio, n + 1)
            end
        end

        for anio_obj in anio:anio_fin
            for ((vec, est, hora), (suma, n)) in acumulados
                precios_mercado_horizonte[(anio_obj, vec, est, hora)] = suma / n
            end
        end
    end
    # Si seguimos sin datos (enero del primer año) → precio plano por defecto
    if isempty(precios_mercado_horizonte)
        vectores   = unique(String.(datos_bundle.demanda.vector))
        estaciones = unique(values(datos_bundle.cal_estaciones))
        for anio_obj in anio:anio_fin, vec in vectores, est in estaciones, hora in 1:24
            precios_mercado_horizonte[(anio_obj, vec, est, hora)] = 50.0
        end
    end

    if ACTIVACION_PREDEF_MILP
        println("    💡 Ejecutando MILP multi-agente (horizonte 5 años).  (para cambiar, ACTIVACION_PREDEF_MILP = false)")
        try
            # 2️⃣  Datos por agente (cuota de mercado & déficits)
            datos_por_agente = Dict{Int,NamedTuple}()
            total_cap        = sum(sum(values(a.capacidades)) for a in values(agentes))

            for (id, ag) in agentes
                cuota = max(0,
                            total_cap > 0 ? sum(values(ag.capacidades)) / total_cap : 0.0)
                demanda_adj = transform(datos_horizonte.demanda,
                                                :valor => (x -> x .* cuota) => :valor)
                # Snapshot de la capacidad instalada (MW) por tecnología y periodo
                hist_cap = Dict{String,Dict{Tuple{Int,Int},Float64}}()
                for (tec, mw) in ag.capacidades
                    # Rellenamos todos los meses del horizonte con la capacidad actual
                    hist_cap[tec] = Dict{Tuple{Int,Int},Float64}(
                        (a, m_) => mw  for a in (anio-1):anio_fin, m_ in 1:12
                    )
                end
                datos_por_agente[id] = (
                    demanda_ajustada     = demanda_adj,
                    deficits_recientes   = deficits_por_vector,
                    cap_existente_agente = hist_cap,
                )
            end

            # 3️⃣  Resolver MILP
            res = OptimizacionSistema.ejecutar_optimizacion_multiagente(
                agentes, datos_por_agente, datos_horizonte,
                Dict(id => tecs for id in keys(agentes)),
                esc, precios_mercado_horizonte;
                tiempo_max = 600
            )

            # 4️⃣  Aplicar sólo cap_nueva del (anio, mes) actual
            aplicar_inversiones_primer_periodo!(
                agentes, res, anio, mes, tecs, datos_bundle.delay_construccion
            )
        catch err
            @warn "MILP falló en $mes/$anio: $err. Usando fallback PlanificacionMes."
            PlanificacionMes.evaluar_inversiones_mes!(
                agentes, tecs, anio, mes, precios_historicos,
                datos_bundle, precio_co2_anio, deficits_por_vector,
                precio_medio_vector
            )
        end
    else
        println("    💡 Ejecutando inversiones con PlanificacionMes (para cambiar, ACTIVACION_PREDEF_MILP = true)")
        PlanificacionMes.evaluar_inversiones_mes!(
            agentes, tecs, anio, mes, precios_historicos,
            datos_bundle, precio_co2_anio, deficits_por_vector,
            precio_medio_vector
        )
    end
end

"""
aplicar_inversiones_primer_periodo!( … )
            
Registra la inversión propuesta para el mes corriente, verificando la  
solvencia del agente y aplicando el delay tecnológico correspondiente.
"""
function aplicar_inversiones_primer_periodo!(
        agentes::Dict{Int,Agentes.EnergyAgent},
        resultado_milp,
        anio_actual::Int,
        mes_actual::Int,
        tecs::Dict{String,Tecnologias.Tecnologia},
        delay_construccion::Dict
)
    for (ag_id, info) in resultado_milp.sub
        if info.opt_status ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
            df_cap = get(info.opt_resultados, :capacidad, DataFrame())
            for row in eachrow(df_cap)
                if row.anio == anio_actual && row.mes == mes_actual && row.cap_nueva > 1e-6
                    if OptimizacionSistema.puede_financiar_inversion(
                                    agentes[ag_id], tecs, row.tecnologia, row.cap_nueva)

                        κ           = get(delay_construccion, row.tecnologia, 0)
                        ready_month = anio_actual*12 + mes_actual + κ
                        InversionesPendientes.registrar_pendiente!(
                                ag_id, ready_month, row.tecnologia, row.cap_nueva
                        )
                        println("      ✅ $(agentes[ag_id].name): $(round(row.cap_nueva,digits=1)) MW " *
                                "en $(row.tecnologia) (delay $(κ) meses)")
                    end
                end
            end
        end
    end
end

# Suma la producción de un día (clave = (ag, tec, hora) o NamedTuple similar)
# a los acumuladores anual / mensual que usa el post-proceso.
function acumular_produccion!(
        prod_anual::Dict{Int, Dict{String, Float64}},
        prod_mensual::Dict{Int, Dict{Int, Dict{String, Float64}}},
        fecha::Date,
        prod_dia::Dict,
)

    anio = year(fecha)
    mes  = month(fecha)
    pa_anio = get!(prod_anual,  anio, Dict{String, Float64}())
    pm_anio = get!(prod_mensual, anio, Dict{Int, Dict{String, Float64}}())
    pm_mes  = get!(pm_anio, mes, Dict{String, Float64}())

    # Recorremos cada sub‐dict (agente → {tec=>MWh}), luego cada (tec, mwh)
    for dia_agente in values(prod_dia)
        for (tec, mwh) in dia_agente
            mwh > 1e-6 || continue                      # ignorar ceros numéricos
            pa_anio[tec] = get(pa_anio, tec, 0.0) + mwh
            pm_mes[tec]  = get(pm_mes,  tec, 0.0) + mwh
        end
    end
end

end # module