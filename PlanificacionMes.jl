module PlanificacionMes
using JuMP
import JuMP: @constraint, @objective
using DataFrames  
using Main.OptimizacionSistema    # para update_mes!, postprocess_mes_capex!, etc.
using Main.Utils                 # para crf()
using Main.Tecnologias: Tecnologia  # importar el tipo Tecnología
using Dates
using ..Agentes
using Main.FinanzasParametros
import Main.OptimizacionSistema: get_cap_existente, get_cap_existente_agente
import JuMP: copy_model
using ..InversionesPendientes  # NUEVO - usar el módulo compartido

using Base.Threads
using HiGHS, Statistics
using DataFrames: groupby

export planificar_mes!, evaluar_inversiones_mes!, evaluar_inversiones_agente!



"""
Evalúa inversiones mensuales para todos los agentes.
Cada agente decide si construir nueva capacidad basándose en rentabilidad esperada.
"""
function evaluar_inversiones_mes!(
    agentes,
    tecnologias,
    anio,
    mes,
    precios_historicos,
    datos_bundle,
    precio_co2_anio,
    deficits_por_vector::Dict{String,Float64},
    precio_medio_vector::Dict{String,Float64}
)
    println("      💡 Evaluando inversiones para mes $mes/$anio")
    
    # Recogemos las parejas (id, agente) en un vector para poder indexar
    pares = collect(agentes)
    Threads.@threads for i in 1:length(pares) 
        ag_id, agente = pares[i]
        evaluar_inversiones_agente!(
            agente, tecnologias, anio, mes,
            precios_historicos, datos_bundle, agentes,
            precio_co2_anio,
            deficits_por_vector,
            precio_medio_vector
        ) # AL LANZARSE AHORA CON THREADS HAY QUE TENER EN CUENTA QUE LOS LOGS SALEN ENTREMEZCLADOS
    end
end

"""
Evalúa inversiones para un agente específico.
"""
function evaluar_inversiones_agente!(
    agente,
    tecnologias,
    anio,
    mes,
    precios_historicos,
    datos_bundle,
    agentes,
    precio_co2_anio,
    deficits_por_vector::Dict{String,Float64}=Dict{String,Float64}(),
    precio_medio_vector::Dict{String,Float64}=Dict{String,Float64}(),
)



    # === CRITERIOS DE ELEGIBILIDAD ESTRICTOS ===
    
    # 1. Liquidez mínima
    #cash_minimo_operacion = 30.0*1e6  # €
    #if agente.cash < cash_minimo_operacion
    #    return
    #end
    
    # 2. Ratio de deuda máximo usando funciones existentes
    ratio_deuda_actual = calcular_ratio_deuda(agente)
    #if ratio_deuda_actual > 3.5  # Aumentado para permitir más inversión
    #    return
    #end

    # 3. Cobertura de cuotas existentes
    #cuotas_anuales = sum(loan.cuota for loan in agente.loans; init=0.0)
    #if cuotas_anuales > agente.revenues * 0.35  # Máximo 35% de ingresos
    #    return
    #end
    
    # 4. Precios de mercado mínimos
    precio_promedio = calcular_precio_promedio_reciente(precios_historicos)
    #if precio_promedio < 25.0  # €/mwh mínimo
    #    return
    #end
    
    # Detectar si hay déficit reciente que justifique inversión de emergencia
    hay_deficit_reciente = false
    tecnologias_deficit = Set{String}()



    
    # Revisar últimos 7 días de precios
    fechas_recientes = sort(collect(keys(precios_historicos)), rev=true)[1:min(7, length(precios_historicos))]
    for fecha in fechas_recientes
        for (hora, vector) in keys(precios_historicos[fecha])
            precio = precios_historicos[fecha][(hora, vector)]
            # Precio >= 500 indica escasez
            if precio >= 500.0
                hay_deficit_reciente = true # ANTIGUAMENTE ERA ASÍ, AHORA SIGUE ESTANDO PORQUE AYUDA A LOS AGENTES A DETECTAR DONDE PUEDEN INVERTIR
                # Buscar tecnologías que producen este vector
                for (tec_name, tec_data) in tecnologias
                    if tec_data.vector == vector
                        push!(tecnologias_deficit, tec_name)
                    end
                end
            end
        end
    end


    
    # 5. Máximo una inversión por trimestre
    #if mes % 3 != 1  # Solo en enero, abril, julio, octubre
    #    return
    #end
    
    # === EVALUACIÓN DE INVERSIONES ===
    
    mejores_opciones = []
    
    println("        🔍 Evaluando inversiones para $(agente.name):")
    println("           - Cash: $(round((agente.cash/1e6), digits=2))M€")
    println("           - Ratio deuda: $(round(ratio_deuda_actual, digits=2))")
    println("           - Hay déficit reciente: $hay_deficit_reciente")
    if hay_deficit_reciente
        println("           - Tecnologías potencialmente deficitarias: $tecnologias_deficit")
    end
    
    # Ver inversiones pendientes de otros agentes
    pend_otros = InversionesPendientes.obtener_pendientes_global()
    delete!(pend_otros, agente.id)
    for (tecnologia, tec_data) in tecnologias
        es_prioritaria = tecnologia in tecnologias_deficit

        if should_consider_technology(agente, tecnologia, tec_data, anio, mes)
            # Para importación con déficit, usar tamaño mayor
            deficit_vector = get(deficits_por_vector, tec_data.vector, 0.0)
        
            if es_prioritaria && occursin("import", lowercase(tecnologia))
                mw_inversion = determine_investment_size_emergency(
                    agente,
                    tecnologia,
                    anio,
                    mes,
                    agentes,
                    deficit_vector
                )
            else
                mw_inversion = determine_investment_size_conservative(agente, tecnologia, anio, mes, agentes)
            end



            demanda_vector = 0.0
            for row in eachrow(datos_bundle.demanda)
                if row.anio == anio && row.vector == tec_data.vector
                    demanda_vector += row.valor
                end
            end

            capex_total = mw_inversion * tec_data.costo_inicial
            
            # Para importación (capex = 0), permitir siempre
            if capex_total < 1e-6 && es_prioritaria
                push!(mejores_opciones, (
                    tecnologia = tecnologia,
                    mw = mw_inversion,
                    capex = 0.0,
                    tir = 999.0  # TIR infinita para capex = 0
                ))
                continue
            end
            
            # Verificar que puede permitirse la inversión
            equity_necesario = capex_total * EQ_SHARE  # Usar constante existente
            #if agente.cash < equity_necesario + cash_minimo_operacion
            #    continue
            #end
            
            if isempty(agente.historial_penalizaciones)
                # primera vez que se usa
                agente.historial_penalizaciones = Dict{Tuple{Int,Int},Float64}()
            end
            # --- penalizaciones medias de los 3 meses previos --------------------------
            refs = meses_anteriores(anio, mes; n = 3)

            vals_pen = [
                get(agente.historial_penalizaciones, par, 0.0) for par in refs
            ]
            
            prom_pen = isempty(vals_pen) ? 0.0 : mean(vals_pen)
            # ---------------------------------------------------------------------------


            # VPN con metodología conservadora
            vpn_esperado = calcular_vpn_conservador_compatible(
                mw_inversion, tecnologia, tec_data, 
                precio_promedio, anio, datos_bundle.perfiles_generacion, precio_co2_anio,
                precio_medio_vector
            )

            # --- NUEVO: ahorro en multas ----------------------------------------------
            pen_vpn = calcular_vpn_ahorro_penalizaciones(
                        mw_inversion,
                        deficit_vector,
                        prom_pen;
                        tasa_descuento = 0.10,
                        horizonte_meses = 60)

            vpn_total = vpn_esperado + pen_vpn
            tir_aproximada = vpn_total / capex_total
            # --------------------------------------------------------------------------



            tir_minima = es_prioritaria ? 0.02 : 0.07
            if hay_deficit_reciente && (tecnologia in tecnologias_deficit)
                tir_minima = 0.02
                if deficit_vector > 0.1 * demanda_vector
                    mw_inversion *= 1.5
                    capex_total = mw_inversion * tec_data.costo_inicial
                    vpn_esperado = calcular_vpn_conservador_compatible(
                        mw_inversion, tecnologia, tec_data,
                        precio_promedio, anio, datos_bundle.perfiles_generacion, precio_co2_anio,
                        precio_medio_vector
                    )
                    # --- NUEVO: ahorro en multas ----------------------------------------------
                    pen_vpn = calcular_vpn_ahorro_penalizaciones(
                        mw_inversion,
                        deficit_vector,
                        prom_pen;               
                        tasa_descuento = 0.10,
                        horizonte_meses = 60)

                    vpn_total = vpn_esperado + pen_vpn
                    tir_aproximada = vpn_total / capex_total
                    # ----------------------
                end
            end


            # --- penalizaciones medias de los 3 meses previos --------------------------
            refs = meses_anteriores(anio, mes; n = 3)

            vals_pen = [
                get(agente.historial_penalizaciones, par, 0.0) for par in refs
            ]
            
            prom_pen = isempty(vals_pen) ? 0.0 : mean(vals_pen)
            # ---------------------------------------------------------------------------
            if prom_pen > 1e6          # > 1 M€/mes de multa media
                println("           ⚠️  $(agente.name) multa media $(round(prom_pen/1e6; digits=2)) M€/mes → se relaja la TIR mínima")
                tir_minima *= 0.5
            end

            # Solo considerar si TIR > 25%
            # Penalizar saturación de tecnologías según inversiones de otros agentes
            cap_futura = sum((inv[3] for v in values(pend_otros) for inv in v if inv[2] == tecnologia); init = 0.0)
            if demanda_vector > 1e-6 && cap_futura > 0
                tir_minima *= (1 + cap_futura / demanda_vector)
            end
            if tir_aproximada > tir_minima
                push!(mejores_opciones, (
                    tecnologia = tecnologia,
                    mw = mw_inversion,
                    capex = capex_total,
                    tir = tir_aproximada
                ))
            end

            
            if mw_inversion > 0.0
                #println("           📊 Evaluando $tecnologia:")
                #println("              - MW inversión: $(round(mw_inversion, digits=1))")
                #println("              - CAPEX: $(round(capex_total, digits=1))€")
                #println("              - TIR: $(round(tir_aproximada*100, digits=1))%")
                #println("              - TIR mínima: $(round(tir_minima*100, digits=1))%")
                # REVISAR, LOG SI HAY PROBLEMAS CON EL REGISTRO DE INVERSIONES, IMPRIMIR RALENTIZA MUCHO EL PROGRAMA
            end


        end

    end
    
    # Agrupar las opciones por vector energético

    # mejores_opciones :: Vector{NamedTuple}
    por_vector = Dict{String,NamedTuple}()

    for opt in mejores_opciones
        v = tecnologias[opt.tecnologia].vector
        if !haskey(por_vector, v) || opt.tir > por_vector[v].tir
            por_vector[v] = opt
        end
    end

    # Ejecutar una inversión por cada vector energético
    for mejor in values(por_vector)
        # Chequea liquidez y ratio de deuda
        if Agentes.registrar_inversion!(agente, mejor.capex; solo_chequear=true)
            Agentes.registrar_inversion!(agente, mejor.capex)

            κ = get(datos_bundle.delay_construccion, mejor.tecnologia, 0)
            InversionesPendientes.registrar_pendiente!(agente.id,
                                                      anio*12 + mes + κ,
                                                      mejor.tecnologia,
                                                      mejor.mw)
            println("      📊 Datos inversión ejecutada: tecnologia = $(mejor.tecnologia), mw = $(mejor.mw), capex = $(mejor.capex), tir = $(mejor.tir)
            Se terminará de instalar en $(κ) meses")

            if !haskey(agente.capacidades, mejor.tecnologia)
                agente.capacidades[mejor.tecnologia] = 0.0
            end
        end
    end
end

# Nueva función auxiliar compatible
function determine_investment_size_conservative(
    agente,
    tecnologia,
    anio,
    mes,
    agentes
)::Float64


    # Prohibir de entrada cualquier inversión en import_electricidad
    if occursin("import_electricidad", lowercase(tecnologia))
        return 0.0
    end

    # (Resto del código idéntico al existente:)
    patrimonio_neto = calcular_patrimonio_neto(agente)
    if agente.cash <= 0 || patrimonio_neto <= 0
        return 0.0
    end

    cap_limit = calcular_cap_limit_trimestral(agente, tecnologia, anio, mes, agentes)
    heuristico = 0.0
    if occursin("import", lowercase(tecnologia))
        if !occursin("import_electricidad", lowercase(tecnologia))
            heuristico = 1e4
        end
    end

    tamanio_estimado = max(heuristico, cap_limit)
    return max(0.0, min(tamanio_estimado, agente.cash * 0.3))
end


function calcular_vpn_conservador_compatible(mw, tecnologia, tec_data, precio_promedio, anio_base, perfiles_generacion_df, precio_co2_anio, precio_medio_vector)
    # Usar parámetros del sistema existente
    tasa_descuento = 0.10  # Más conservador
    horizonte_anios = 15
    
    factor_capacidad = average_capacity_factor(tecnologia, perfiles_generacion_df)  # Función existente, REVISAR E INCORPORAR PERFILES
    mwh_anuales = mw * 8760 * factor_capacidad
    
    vpn = 0.0
    for anio in 1:horizonte_anios
        # Degradación de precios por competencia
        precio_anio = get(precio_medio_vector, tecnologia, precio_promedio) *
                      (0.98)^(anio - 1)
        ingresos_anuales = mwh_anuales * precio_anio 
        
        # Costos usando estructura existente
        om_fijo_anual = mw * tec_data.costo_om
        om_variable_anual = mwh_anuales * tec_data.costo_om_variable
        
        # ─── Coste del vector consumido ────────────────────────────────────────
        # Solo aplica cuando `vector_consumido` no está vacío y el ratio es > 0.
        costo_energia = 0.0
        if !isempty(tec_data.vector_consumido) && tec_data.ratio_energia > 0
            costo_energia = get(precio_medio_vector,
                                tec_data.vector_consumido,   # devuelve 0.0 si no existe
                                0.0) *
                            mwh_anuales / tec_data.ratio_energia
        end
        
        flujo_anual = ingresos_anuales - om_fijo_anual - om_variable_anual -
                      tec_data.emisiones_unitarias * precio_co2_anio * mwh_anuales -
                      costo_energia
        if tec_data.emisiones_unitarias == 0
            # Asumir subvención promedio de 40 €/MWh para renovables  
            subvencion_esperada = mwh_anuales * 40.0
            flujo_anual += subvencion_esperada
        end

        vpn += flujo_anual / (1 + tasa_descuento)^anio
    end
    
    return vpn
end

"""
    calcular_vpn_ahorro_penalizaciones(
        mw_inversion, deficit_vector, prom_pen_mensual;
        tasa_descuento = 0.10, horizonte_meses = 60)

Devuelve el VPN del ahorro en multas que se obtendría al cubrir
parte del déficit del *vector* durante los próximos `horizonte_meses`.
"""
function calcular_vpn_ahorro_penalizaciones(
        mw_inversion::Float64,
        deficit_vector::Float64,
        prom_pen_mensual::Float64;
        tasa_descuento::Float64 = 0.10,
        horizonte_meses::Int   = 60)

    # fracción del déficit que cubro con la nueva capacidad
    fraccion_cubierta = deficit_vector > 0 ?
                        min(1.0, mw_inversion / deficit_vector) : 0.0

    ahorro_mensual = prom_pen_mensual * fraccion_cubierta
    vpn = 0.0
    for mes in 1:horizonte_meses
        vpn += ahorro_mensual / (1 + tasa_descuento)^(mes/12)
    end
    return vpn
end



"""
    average_capacity_factor(tecnologia, perfiles_df)

Calcula el factor de capacidad promedio para `tecnologia` usando
`perfiles_df` (formato [:tecnologia, :estacion, :hora, :factor]).
Devuelve `1.0` si la tecnología no está presente.
"""
function average_capacity_factor(tecnologia::String, perfiles_df)
    if isempty(perfiles_df) || :tecnologia ∉ names(perfiles_df)
        return 1.0
    end
    sub = perfiles_df[perfiles_df.tecnologia .== tecnologia, :factor]
    return isempty(sub) ? 1.0 : mean(sub)
end



function calcular_precio_promedio_reciente(precios_historicos)
    if isempty(precios_historicos)
        return 50.0  # Precio por defecto
    end
    
    # Promedio de los últimos 30 días
    fechas_recientes = sort(collect(keys(precios_historicos)), rev=true)[1:min(30, length(precios_historicos))]
    
    total_precios = 0.0
    contador = 0
    
    for fecha in fechas_recientes
        for (hora, precio) in precios_historicos[fecha]
            total_precios += precio
            contador += 1
        end
    end
    
    return contador > 0 ? total_precios / contador : 50.0
end

function should_consider_technology(agente, tecnologia, tec_data, anio, mes)::Bool
    # Filtros básicos para decidir si evaluar una tecnología
    
    # ───────────────── Nuevo filtro de capacidad global ─────────────────
    # Si con esta inversión la capacidad total instalada del sistema
    # superase los 2 000 000 MW, se descarta.
    cap_total_sistema = get_cap_existente(tecnologia, anio, mes)
    if (cap_total_sistema > 200_000.0 || occursin("import", lowercase(tecnologia)))    # GET_CAP_EXISTENTE NO FUNCIONA BIEN, ASÍ QUE CORTO POR LO SANO
        #println("        ⚠️  Inversión descartada para $(tecnologia): "
        #                * "$(round(cap_total_sistema, digits = 0)) MW "
        #                * "superarían el límite de 2 000 000 MW.")
        return false
    end


    # No invertir en tecnologías que ya tiene en exceso
    #capacidad_actual = get(agente.capacidades, tecnologia, 0.0)
    #if capacidad_actual > 1000000.0  # Más de 1000000 mw
    #    return false
    #end
    
    # No invertir en tecnologías muy contaminantes /si es un escenario verde
    if tec_data.emisiones_unitarias > 500.0  # Más de 500 toneladas co2/mwh
        return false
    end
    
    # No invertir en importaciones eléctricas
    if occursin("import", lowercase(tecnologia))
        if occursin("import_electricidad", lowercase(tecnologia))
            return false    #NO SE PUEDE AUMENTAR LA IMPORTACION DE ELECTRICIDAD, ES ASUNTO ESTATAL
        else
            return true
        end
    end
    
    return true
end


function determine_investment_size_emergency(
    agente::Any,         # EnergyAgent
    tecnologia::String,
    anio::Int,
    mes::Int,
    agentes,
    deficit_vector::Float64
)::Float64
    # 0) Prohibir de entrada cualquier inversión en import_electricidad
    if occursin("import_electricidad", lowercase(tecnologia))
        return 0.0
    end

    # 1) Validaciones financieras (cash > 0, patrimonio > 0)
    patrimonio_neto = calcular_patrimonio_neto(agente)
    if agente.cash <= 0 || patrimonio_neto <= 0
        return 0.0
    end

    # 2) Calcular el límite estricto para este mes:
    cap_limit = calcular_cap_limit_trimestral(agente, tecnologia, anio, mes, agentes) * 1.5

    # 3) HEURÍSTICO PARA IMPORTACIONES (solo para importaciones que NO sean import_electricidad)
    heuristico = 0.0
    if occursin("import", lowercase(tecnologia))
        if !occursin("import_electricidad", lowercase(tecnologia))
            heuristico = max(1e6, deficit_vector * 1.2)
        end
    end

    # 4) Limitar el tamaño final al mínimo entre “cap_limit” y “heurístico”,
    #    y también condicionarlo a su liquidez (80% del cash).
    #    (Ya hemos bloqueado import_electricidad, así que aquí no hace falta re-chequearlo.)
    tamanio_estimado = max(heuristico, cap_limit)
    return max(0.0, min(tamanio_estimado, agente.cash * 0.8))
end

"""
dado un agente, una tecnología y un (año, mes), devuelve el
límite de MW que puede invertir ese agente en ese mes,
según la regla “4% de lo instalado o 1000 de mínimo”.

Estaba en PlanificacionMes.jl
SE PUEDE RESCATAR LO DEL % DE INSTALACIÓN YA INSTALADA PARA SIMULAR LOS CRECIMIENTOS POCO A POCO DE LA VIDA REAL
"""
function calcular_cap_limit_trimestral(agente, tecnologia, anio, mes_real, agentes)

    # 1) Construimos el Periodo correspondiente al mes anterior.
    #    Suponemos que Periodo se construye como Periodo(año, mes).
    mes_prev = mes_real - 1
    anio_prev = anio
    if mes_prev == 0
        mes_prev = 12
        anio_prev -= 1
    end
    periodo_prev = (anio_prev, mes_prev)

    cap_inst = get_cap_existente_agente(agente.id, tecnologia, periodo_prev, agentes)
    # Si quisiéramos evitar que un agente instale 100 MW cuando no tenga nada instalado, 
    # podríamos cambiar a: if cap_inst == 0 return 0 else max(0.04*cap_inst, 100) end.
    return max(0.08 * cap_inst, 150.0)

end


"""
Devuelve tuplas (año, mes) de los `n` meses inmediatamente anteriores
al par (`año`,`mes_actual`), sin incluir el mes actual.
"""
function meses_anteriores(anio::Int, mes_actual::Int; n::Int = 3)
    refs = Tuple{Int,Int}[]
    a, m = anio, mes_actual
    for _ in 1:n
        m -= 1
        if m == 0          # pasamos de enero a diciembre del año anterior
            a -= 1
            m = 12
        end
        push!(refs, (a, m))
    end
    return refs
end


end # module
