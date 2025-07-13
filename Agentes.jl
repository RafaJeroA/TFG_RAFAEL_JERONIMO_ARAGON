module Agentes

using CSV, DataFrames
using Dates
using ..FinanzasParametros
using ..Utils

cuota_constante(P; i, n) = i ≈ 0 ? P / n : P * (i * (1 + i)^n) / ((1 + i)^n - 1)

"""
Módulo Agentes: definición de EnergyAgent y funciones ABM.

Ejemplo:
    # Supongamos un agente `ag` con capacidades existentes
    req = calcular_expansion_deseada(ag, ag.capacidades)
"""



# Intento de carga de Agents.jl para obtener implementaciones reales
try
    import Agents: AbstractAgent, ABM, run!
    @info "Agents.jl cargado correctamente"
catch e
    @warn "Agents.jl no disponible: $e. Usando stub ABM."
    # Stubs ABM por defecto si Agents.jl no está disponible
    struct MissingABMModel end
    abstract type AbstractAgent end
    function ABM(::Type; kwargs...); MissingABMModel(); end
    function run!(::MissingABMModel, steps); end
end


const _SEQ_PRESTAMO = Ref(0)
next_loan_id() = (_SEQ_PRESTAMO[] += 1)


mutable struct Loan
    id_prestamo::Int
    id_agente_deudor::Int
    principal::Float64
    saldo::Float64
    plazo::Float64     # anios restantes
    interes::Float64   # tipo anual
    cuota::Float64     # cuota anual
    
    # Constructor interno para validar tipos
    function Loan(id_prestamo::Int, id_agente_deudor::Int, principal::Real, saldo::Real, plazo::Real, interes::Real, cuota::Real)
        new(Int(id_prestamo), Int(id_agente_deudor), Float64(principal), Float64(saldo), Float64(plazo), Float64(interes), Float64(cuota))
    end
end
#println("DEBUG: Definición de Agentes.Loan cargada con objectid: ", objectid(Loan)) # <--- AÑADE ESTO


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Definición del tipo de agente                                           ║
# ╚══════════════════════════════════════════════════════════════════════════╝
"""
Un agente energético con rol (`:producer` o `:consumer`), capacidad/demanda y
precio.  
**Nota**: `id` lo asigna automáticamente `Agents.jl`; no lo seteamos a mano.
"""
mutable struct EnergyAgent <: AbstractAgent
    id::Int
    name::String
    capacidades::Dict{String,Float64}  # mw por tecnología
    var_costs::Dict{String,Float64}   # €/mwh por tecnología
    emissions::Dict{String,Float64}   # tco2/mwh por tecnología
    cash::Float64
    debt::Float64
    revenues::Float64
    profit::Float64
    cost::Float64
    price::Float64      # precio ofertado/pagado
    equity::Float64              # €
    opex_fijo::Float64           # €/año (suma cost_OM * mw)
    capex_acum::Float64          # €
    EBITDA::Float64              # €/año   
    loans::Vector{Loan}
    produccion_anual::Dict{String,Float64}
    produccion_mensual::Dict{Int,Dict{Int,Dict{String,Float64}}}
    deficits_mensual::Dict{Int,Dict{Int,Dict{String,Float64}}}  
    produccion_diaria::Dict{Date,Dict{String,Float64}}
    profit_mes::Float64              # /mes - beneficio del mes actual
    emissions_mes::Float64           # tco2/mes - emisiones del mes actual
    historial_penalizaciones::Dict{Tuple{Int,Int},Float64}
end






const CASH_RATIO = 0.50   # 50 % del beneficio va a caja
const DEBT_RATIO = 0.50   # 50 % del coste variable se financia como deuda



#"""
#    calcular_expansion_deseada(ag::EnergyAgent, cap_cs::Dict{String,Float64}, max_pct::Float64 = MAX_EXPANSION_PCT)
#Devuelve Dict{String,Float64} con mw extra por tecnología.
#"""
#function calcular_expansion_deseada(ag::EnergyAgent, cap_cs::Dict{String,Float64}, max_pct::Float64 = MAX_EXPANSION_PCT)
#    req = Dict{String,Float64}()
#    for (tecnologia, cap) in cap_cs
#        req[tecnologia] = min(max_pct * cap, cap)
#    end
#    return req
#end

# Exported symbols
export EnergyAgent, inicializar_model,
       cargar_agentes_generacion, cargar_param_agentes, guardar_finanzas,
       cargar_finanzas_iniciales!,
       EBITDA, calcular_patrimonio_neto, calcular_ratio_deuda, calcular_expansion_deseada, 
       registrar_inversion!, es_inversion_rentable, puede_servir_deuda,
       determine_investment_size_conservative





function EBITDA(ag::EnergyAgent)
    return ag.revenues - ag.cost - ag.opex_fijo
end

"""
cargar_agentes_generacion(file)
Lee CSV de agentes de generación y retorna Dict{agent_id=>EnergyAgent}.
"""
function cargar_agentes_generacion(file::String)
    # Leer capacidades de agentes
    df_agents = CSV.read(file, DataFrame)
    # Leer datos de tecnologías para costes variables y emisiones
    tec_df = CSV.read(joinpath(@__DIR__, "Datos", "tecnologias.csv"), DataFrame)
    # Construir mapas de coste variable y emisiones
    var_cost_map = Dict{String,Float64}()
    emis_map = Dict{String,Float64}()
    for row in eachrow(tec_df)
        t_original = String(row.tecnologia)  # Nombre original del CSV
        t_normalizado = normalizar_nombre_tecnologia(t_original)
        costo_var = Float64(coalesce(get(row, :costo_om_variable, 0.0), 0.0))
        emis_val = Float64(coalesce(get(row, :emisiones_unitarias, 0.0), 0.0))
        
        # Guardar con ambos nombres para compatibilidad
        var_cost_map[t_original] = costo_var
        var_cost_map[t_normalizado] = costo_var
        emis_map[t_original] = emis_val
        emis_map[t_normalizado] = emis_val
    end
    
    agents = Dict{Int,EnergyAgent}()
    for r in eachrow(df_agents)
        aid = Int(r.agent_id)
        if !haskey(agents, aid)
            agents[aid] = EnergyAgent(aid, String(r.empresa),
                               Dict(), Dict(), Dict(),
                               1e10,   # cash
                               0.0,   # debt
                               0.0,   # revenues
                               0.0,   # profit
                               0.0,   # cost
                               0.0,   # price
                               0.0,   # equity
                               0.0,   # opex_fijo
                               0.0,   # capex_acum
                               0.0,   # EBITDA
                               Vector{Loan}(),
                               Dict{String,Float64}(),
                               Dict{Int,Dict{String,Float64}}(),
                               Dict{Int,Dict{String,Float64}}(),
                               Dict{Date,Dict{String,Float64}}(),
                               0.0,   # profit_mes
                               0.0,   # emissions_mes
                               Dict{Tuple{Int,Int},Float64}()
                               )
        end
        ag = agents[aid]
        tecnologia = normalizar_nombre_tecnologia(String(r.tecnologia))
        ag.capacidades[tecnologia] = Float64(r.capacidad_mw)
        ag.var_costs[tecnologia]   = get(var_cost_map, tecnologia, 0.0)
        ag.emissions[tecnologia]   = get(emis_map,    tecnologia, 0.0)   # ← usa emissions
    end

    # === VERIFICACIÓN Y SINCRONIZACIÓN FINAL ===
    println("🔧 Verificando sincronización de costos para todos los agentes...")
    for (ag_id, agente) in agents
        tecnologias_faltantes = String[]
        
        for tecnologia in keys(agente.capacidades)
            if !haskey(agente.var_costs, tecnologia)
                push!(tecnologias_faltantes, tecnologia)
                # Buscar en var_cost_map (que viene de tecnologias.csv)
                var_cost_valor = get(var_cost_map, tecnologia, 0.0)
                agente.var_costs[tecnologia] = var_cost_valor
                
                # Hacer lo mismo para emisiones
                emis_valor = get(emis_map, tecnologia, 0.0)
                agente.emissions[tecnologia] = emis_valor
            end
        end
        
        if !isempty(tecnologias_faltantes)
            println("  ⚠️  Agente $(agente.name): Sincronizadas $(length(tecnologias_faltantes)) tecnologías: $tecnologias_faltantes")
        end
    end


    for ag in values(agents)
        ag.capacidades = Dict(
            Utils.normalizar_nombre_tecnologia(k) => v
            for (k,v) in ag.capacidades
        )
        ag.var_costs = Dict(
            Utils.normalizar_nombre_tecnologia(k) => v
            for (k,v) in ag.var_costs
        )
        ag.emissions = Dict(                              # ← emissions
            Utils.normalizar_nombre_tecnologia(k) => v
            for (k,v) in ag.emissions
        )
    end
    


    return agents
end



# === Finanzas Iniciales y Resultados ===


"""
guardar_finanzas(agents, scenario)
Escribe CSV con finanzas tras la optimización para cada agente.
"""
function guardar_finanzas(agents::Dict{Int,EnergyAgent}, scenario::String)
    timestamp = Dates.format(now(), "yyyy_mmdd_HHMMSS")
    rows = [ (ag.id, ag.name, ag.cash, ag.debt, ag.revenues, ag.profit) for ag in values(agents) ]
    df = DataFrame(agent_id=getindex.(rows,1),
                   empresa=getindex.(rows,2),
                   cash_eur=getindex.(rows,3),
                   debt_eur=getindex.(rows,4),
                   revenues_eur=getindex.(rows,5),
                   profit_eur=getindex.(rows,6))
    path = joinpath(@__DIR__, "Datos", "agentes_finanzas_$(scenario)_$(timestamp).csv")
    CSV.write(path, df)
    return path
end




function actualizar_finanzas!(ag::EnergyAgent,
                              ingresos::Float64,
                              costo_var::Float64,
                              capex_nuevo::Float64,
                              vida_util::Int)

    # 1. Registrar nuevo préstamo/equity si hay inversión
    if capex_nuevo > 1e-3
        equity = EQ_SHARE * capex_nuevo
        deuda  = capex_nuevo - equity
        cuota  = cuota_constante(deuda; i = INT_RATE, n = PLAZO_PRESTAMO)
        nuevo_loan_af = Loan(
            Int(next_loan_id()),
            Int(ag.id),
            Float64(deuda),
            Float64(deuda),
            Float64(PLAZO_PRESTAMO),
            Float64(INT_RATE),
            Float64(cuota)
        )
        push!(ag.loans, nuevo_loan_af)
    end
    # 2. Costes fijos O&M anuales
    opex_fijo = ag.opex_fijo            # ya en €
    # 3. Depreciación lineal
    deprec = capex_nuevo / vida_util
    # 4. Intereses del periodo
    intereses = 0.0
    for l in ag.loans
        int = l.saldo * l.interes
        amort = l.cuota - int
        l.saldo -= amort
        l.plazo -= 1 #años restantes
        intereses += int
        ag.debt   -= amort
        ag.cash   -= l.cuota
    end
    ag.loans = filter(l -> l.saldo > 1e-3, ag.loans)

    # 5. Cuenta de resultados
    EBIT  = ingresos - costo_var - opex_fijo - deprec
    EBT   = EBIT - intereses
    imp   = max(0, EBT) * TAX_RATE
    ag.profit = EBT - imp
    CFO   = EBIT - imp + deprec        # aproximación
    ag.cash += CFO

    # 6. Actualizar ratios
    ag.revenues = ingresos
end

# -------------------- Finanzas iniciales --------------------

export cargar_finanzas_iniciales!

"""
    cargar_finanzas_iniciales!(agents;
                               ruta = joinpath(@__DIR__, "..", "Datos",
                                               "agentes_finanzas.csv"))

Lee el CSV de finanzas iniciales y:
  * actualiza cash, debt, revenues, cost, profit, equity
  * crea un `Loan` (con cuota constante) si la deuda > 0
"""
function cargar_finanzas_iniciales!(agents::Dict{Int,EnergyAgent};
                                    ruta = joinpath(@__DIR__, "..", "Datos",
                                                    "agentes_finanzas.csv"))
    df = CSV.read(ruta, DataFrame)

    for r in eachrow(df)
        ag = agents[r.agent_id]

        # --- saldos base ---
        ag.cash     = r.cash_eur
        ag.debt     = r.debt_eur
        ag.revenues = r.revenues_eur
        ag.cost     = r.cost_eur
        ag.profit   = r.profit_eur

        # --- equity (aprox.) ---
        ag.equity = ag.cash - ag.debt

        # --- préstamo si procede ---
        # Antes de crear préstamos, sincronizar contador
        if r.debt_eur > 1
            loan_id = next_loan_id()  # ID secuencial consistente
            cuota = cuota_constante(r.debt_eur; i = INT_RATE, n = PLAZO_PRESTAMO)
            nuevo_loan_cfi = Loan(
                Int(loan_id),
                Int(ag.id),
                Float64(r.debt_eur),
                Float64(r.debt_eur),
                Float64(PLAZO_PRESTAMO),
                Float64(INT_RATE),
                Float64(cuota)
            )
            push!(ag.loans, nuevo_loan_cfi)
        end
    end
    return nothing
end

function es_inversion_rentable(agente, tecnologia, mw_nuevos, precio_esperado, tec_data)
    capex = mw_nuevos * tec_data.costo_inicial 
    
    # Ingresos anuales esperados
    factor_capacidad = 1  # CAMBIAR EN FUTURO, METER UNA NUEVA VARIABLE EN EL STRUCT DE TECNOLOGIAS QUE SEA EL PERFIL HORARIO DE CADA UNA PARA MANEJARLO MÁS FÁCIL.
    mwh_anuales = mw_nuevos * 8760 * factor_capacidad
    ingresos_anuales = mwh_anuales * precio_esperado  # €
    
    # Costos anuales
    om_fijo = mw_nuevos * tec_data.costo_om
    om_variable = mwh_anuales * tec_data.costo_om_variable
    
    # VPN simple a 10 años
    flujo_anual = ingresos_anuales - om_fijo - om_variable
    vpn = sum(flujo_anual / (1.08)^anio for anio in 1:10) - capex
    
    return vpn > 0  # Invertir solo si VPN > 0
end


function puede_servir_deuda(agente, nueva_deuda)
    nueva_cuota = cuota_constante(nueva_deuda; i=INT_RATE, n=PLAZO_PRESTAMO)
    cuotas_totales = sum(l.cuota for l in agente.loans) + nueva_cuota
    
    # Usar el nuevo criterio de patrimonio neto en lugar de EBITDA
    ratio_deuda_futuro = calcular_ratio_deuda(agente) + (nueva_deuda / calcular_patrimonio_neto(agente))
    
    return ratio_deuda_futuro < 3.0  # Máximo 300% de ratio de deuda
end

function registrar_inversion!(agente, capex::Float64; solo_chequear::Bool=false)


    if solo_chequear
        println("    📋 Evaluando inversión de $(round(capex, digits=2)) € para $(agente.name)")
    end


    if capex < 1
        return false  # No hay inversión real
    end
    
    equity = EQ_SHARE * capex
    deuda = capex - equity
    
    # ✅ NUEVO CRITERIO: Ratio de deuda basado en patrimonio neto
    # Patrimonio neto = Capacidad total instalada × 7,889,400€/mw (60€/mwh × 24h × 365.25 días × 15 años)
    capacidad_total_mw = sum((v for (k,v) in agente.capacidades); init = 0.0)
    capacidad_propia_mw = sum((v for (k,v) in agente.capacidades if !startswith(k, "import_")); init = 0.0)

    patrimonio_neto = capacidad_propia_mw * 7.8894 * 1e6 + (capacidad_total_mw - capacidad_propia_mw) * 1000
    
    # Ratio de deuda = Total pasivos / Patrimonio neto
    deuda_total_futura = agente.debt + deuda
    ratio_deuda = patrimonio_neto > 0 ? deuda_total_futura / patrimonio_neto : Inf
    
    # ✅ CRITERIO PRINCIPAL: Límite de ratio de deuda
    #max_ratio_deuda = 3.0  # Máximo 300% de deuda sobre patrimonio neto
    #if ratio_deuda > max_ratio_deuda
    #    @warn "$(agente.name): Ratio de deuda demasiado alto ($(round(ratio_deuda, digits=2)) > $max_ratio_deuda). Inversión rechazada."
    #    return
    #end
    # ✅ CRITERIOS MEJORADOS DE INVERSIÓN
    #max_ratio_deuda = 3.5  # Reducido de 3.0 a 2.5
    #cash_minimo_post_inversion = 25.0*1e3  # € mínimo tras inversión
    # AHORA MISMO LOS CRITERIOS QUE HABÍA AQUÍ ESTÁN COMENTADOS PORQUE YA SE APLICAN DIRECTAMENTE EN EL MILP


    if solo_chequear
        return true
    end
    # ✅ EJECUTAR INVERSIÓN
    agente.cash -= equity
    agente.debt += deuda
    agente.capex_acum += capex
    
    # Crear préstamo
    if deuda > 1
        cuota_calculada = cuota_constante(deuda; i=INT_RATE, n=PLAZO_PRESTAMO)
        nuevo_loan = Loan(
            Int(next_loan_id()),
            Int(agente.id),
            Float64(deuda),
            Float64(deuda),
            Float64(PLAZO_PRESTAMO),
            Float64(INT_RATE),
            Float64(cuota_calculada)
        )


        #println("DEBUG registrar_inversion!: El tipo Loan esperado por agente.loans es: ", eltype(agente.loans), " con objectid: ", objectid(eltype(agente.loans)))
        #println("DEBUG registrar_inversion!: El tipo Loan de nuevo_loan es:       ", typeof(nuevo_loan), " con objectid: ", objectid(typeof(nuevo_loan)))

        if objectid(eltype(agente.loans)) !== objectid(typeof(nuevo_loan))
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("ERROR DETECTADO: Los objectid de Loan SON DIFERENTES. Esto confirma la redefinición del tipo.")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        else
            #println("INFO: Los objectid de Loan PARECEN SER IGUALES.")
        end




        push!(agente.loans, nuevo_loan)
        println("      ✅ Inversión ejecutada: $(round(capex, digits=2)) € (Equity: $(round(equity, digits=2)), Deuda: $(round(deuda, digits=2)))")
        println("      📊 Ratio de deuda resultante: $(round(ratio_deuda, digits=2))")
    end
end

function registrar_operacion!(ag, ingresos, cost_var, emis)
    ag.revenues += ingresos
    ag.cost     += cost_var
    ag.profit   += ingresos - cost_var
    ag.emissions_mes += emis
    ag.profit_mes += (ingresos - cost_var)
    ag.cash += ingresos - cost_var
end



"""
Actualiza finanzas mensuales del agente.
"""
function actualizar_finanzas_mes!(agente, anio::Int, mes::Int)
    # Resetear contadores mensuales al inicio
    agente.profit_mes = 0.0
    agente.emissions_mes = 0.0
    
    # O&M fijo mensual
    opex_mensual = agente.opex_fijo / 12
    agente.cash -= opex_mensual
    agente.cost += opex_mensual
    agente.profit_mes -= opex_mensual
    
    # Agregar producción mensual al acumulado anual
    if !hasproperty(agente, :produccion_anual)
        agente.produccion_anual = Dict{String,Float64}()
    end
    
    if hasproperty(agente, :produccion_diaria)
        for (fecha, prod_dia) in agente.produccion_diaria
            if year(fecha) == anio && month(fecha) == mes
                for (tecnologia, mwh) in prod_dia
                    agente.produccion_anual[tecnologia] = get(agente.produccion_anual, tecnologia, 0.0) + mwh  # EN mwh
                end
            end
        end
    end
end

"""
Paga cuotas mensuales de préstamos.
"""
function pagar_cuotas_mensuales!(agentes)
    for ag in values(agentes)
        total_cuotas_mensuales = 0.0
        
        # Calcular total de cuotas mensuales
        for loan in ag.loans
            total_cuotas_mensuales += loan.cuota / 12
        end
        
        if total_cuotas_mensuales < 1
            continue  # No hay cuotas que pagar
        end
        
        if ag.cash >= total_cuotas_mensuales
            # Puede pagar todas las cuotas normalmente
            for loan in ag.loans
                cuota_mensual = loan.cuota / 12
                interes_mensual = loan.saldo * (loan.interes / 12)
                principal_mensual = cuota_mensual - interes_mensual
                
                ag.cash -= cuota_mensual
                loan.saldo -= principal_mensual
                loan.plazo -= 1/12
                ag.debt -= principal_mensual
            end
        elseif ag.cash >= total_cuotas_mensuales * 0.3  # Puede pagar al menos 30%
            # Reestructuración automática
            pago_disponible = ag.cash * 0.8  # Deja 20% de liquidez
            proporcion_pago = min(1.0, pago_disponible / total_cuotas_mensuales)
            
            for loan in ag.loans
                cuota_mensual = loan.cuota / 12
                pago_real = cuota_mensual * proporcion_pago
                
                # Extender plazo proporcionalmente
                extension_meses = (1 - proporcion_pago) * 12
                loan.plazo += extension_meses / 12
                
                # Reducir cuota futura
                loan.cuota *= 0.98  # 2% de reducción
                
                ag.cash -= pago_real
            end
            
            @warn "$(ag.name) reestructura préstamos (pago: $(round(proporcion_pago*100, digits=1))%)"
        else
            # Quita parcial de deuda (bankruptcy protection)
            quita_porcentaje = 0.15  # 15% de quita
            
            for loan in ag.loans
                quita_loan = loan.saldo * quita_porcentaje
                loan.saldo -= quita_loan
                loan.principal -= quita_loan
                ag.debt -= quita_loan
                
                # Recalcular cuota con nuevo principal
                if loan.plazo > 0
                    loan.cuota = cuota_constante(loan.saldo; i=loan.interes, n=loan.plazo)
                end
            end
            
            @warn "$(ag.name) aplica quita de deuda del $(quita_porcentaje*100)%"
        end
        
        # Eliminar préstamos completamente pagados
        filter!(loan -> loan.saldo > 1 && loan.plazo > 1, ag.loans)
    end
end

"""
Calcula el patrimonio neto del agente basado en su capacidad instalada.
Valor estándar: 60€/mwh × 24h × 365.25 días × 15 años = 7,889,400€/mw
"""
function calcular_patrimonio_neto(agente::EnergyAgent)::Float64
    # Sólo sumamos capacidades propias, no las importaciones, ya que esto cuenta como un contrato de energía, no energía que crean estos.
    capacidad_total_mw = sum((v for (k,v) in agente.capacidades); init = 0.0)
    capacidad_propia_mw = sum((v for (k,v) in agente.capacidades if !startswith(k, "import_")); init = 0.0)

    patrimonio_neto = capacidad_propia_mw * 7.8894 * 1e6 + (capacidad_total_mw - capacidad_propia_mw) * 1000
    return patrimonio_neto
end

"""
Calcula el ratio de deuda del agente (Deuda total / Patrimonio neto).
"""
function calcular_ratio_deuda(agente::EnergyAgent)::Float64
    patrimonio = calcular_patrimonio_neto(agente)
    return patrimonio > 0 ? agente.debt / patrimonio : Inf
end

# En Agentes.jl, después de la definición del struct EnergyAgent, añadir:

"""
Calcula el riesgo de penalización para un agente basado en su historial
"""
function calcular_riesgo_penalizacion(agente::EnergyAgent, anio::Int, mes::Int)::Float64
    if isempty(agente.historial_penalizaciones)
        return 0.0
    end
    
    # Considerar últimos 3 meses
    pen_total = 0.0
    meses_considerados = 0
    
    for delta_mes in 0:2
        mes_check = mes - delta_mes
        anio_check = anio
        if mes_check <= 0
            mes_check += 12
            anio_check -= 1
        end
        
        pen = get(agente.historial_penalizaciones, (anio_check, mes_check), 0.0)
        if pen > 0
            pen_total += pen
            meses_considerados += 1
        end
    end
    
    return meses_considerados > 0 ? pen_total / meses_considerados : 0.0
end

export calcular_riesgo_penalizacion


# ──────────────────────────────────────────────────────────────────────────
# NUEVO: registrar_operacion_diaria!
# Registra de forma consistente los resultados de un despacho horario.
# ──────────────────────────────────────────────────────────────────────────

"""
    registrar_operacion_diaria!(agente,
                                mwh_producidos,
                                precio_venta,
                                costo_om_variable,
                                costo_emisiones,
                                subvencion_unitaria = 0.0)
Actualiza todas las magnitudes financieras de *agente* derivadas de la venta diaria:

* **Ingresos**  = *mwh_producidos* × ( `precio_venta` + `subvencion_unitaria` )
* **Costes**    = *mwh_producidos* × ( `costo_om_variable` + `costo_emisiones` )
* **Beneficio** = Ingresos – Costes  →  se acumula en `profit` y `profit_mes`
* **Caja**      += Beneficio
"""
function registrar_operacion_diaria!(
    ag::EnergyAgent,
    mwh::Float64,
    precio_venta::Float64,
    costo_om::Float64,
    costo_emis::Float64,
    subv::Float64 = 0.0,
)
    ingresos = mwh * (precio_venta + subv)
    coste    = mwh * (costo_om + costo_emis)
    emis_t   = mwh * (costo_emis > 0 ? (costo_emis / precio_venta) : 0.0)  # aprox.

    # Reutilizamos la función existente para mantener la coherencia interna
    registrar_operacion!(ag, ingresos, coste, emis_t)

    ag.cash      += ingresos - coste          # caja a día de hoy
    ag.profit_mes += ingresos - coste         # beneficio del mes
end

export registrar_operacion_diaria!


end # module Agentes
