module EquilibrioMultivector

using JuMP, Complementarity, NLsolve, HiGHS
#include("Datos.jl")
using ..Datos, ..Tecnologias
using Logging
using SparseArrays
using MathOptInterface
using SparseDiffTools
using Statistics

using Dates: day, month

using Base.Threads

#global_logger(ConsoleLogger(stderr, Logging.Debug))
#const FINE_DEFICIT  = 3000.0        # €/MWh – multa por déficit
const PENAL_DEF = PARAMS[:precio_base_deficit]

#const DEFICIT_TOL   = 10             # MWh  – tolerancia de déficit

const MOI = MathOptInterface

const MIN_SYNC_SHARE = 0.20 # DADA POR ORGANIZACIONES NACIONALES #REF ADJUNTAR REF


"""
subastar_multivector!(
        fecha,
        agentes, tecnologias,
        demand_ext::Dict,           # Dict{Tuple{Symbol,Int},Float64}
        cap::Dict,                  # Dict{Tuple{Int,Symbol,Int},Float64}
        costo_marg::Dict,           # Dict{Tuple{Int,Symbol},Float64}
        x0 = nothing)               # Semilla para el solver, se utiliza el día de antes, al ser, salvo cambio estacional y meteorológico, muy parecidas las condiciones diarias

Devuelve (q::Dict, P_final::Dict).  Todos los vectores entran en un único
equilibrio Cournot con penalización escalonada por déficit.
"""
function subastar_multivector!(fecha,
                               agentes, tecnologias,
                               demand_ext, cap, costo_marg; x0 = nothing)

    # ----------------------------------------------------------
    # 1. Listas de índices
    # ----------------------------------------------------------
    AGENTES  = collect(keys(agentes))
    raw_vprod   = getfield.(values(tecnologias), :vector)
    VPROD   = unique(Symbol.(raw_vprod))
    raw_vconsum = getfield.(values(tecnologias), :vector_consumido)
    VCONSUM = unique(Symbol.(filter(x -> x != "", raw_vconsum)))
    VECTORES = union(VPROD, VCONSUM)                  # todos los que aparecen

    # ----------------------------------------------------------
    # 1-bis  Parámetros de la inversa de demanda  P(Q)=a_v−b_v·Q
    #        (Dict{Symbol,Float64} en Datos.PARAMS)
    # ----------------------------------------------------------
    a_v = Datos.PARAMS[:dem_inv_a]     # intercepto
    b_v = Datos.PARAMS[:dem_inv_b]     # pendiente   (>0)
    HORAS    = 0:23

    # ----------------------------------------------------------
    # 1-ter  Pendiente de penalización (muy empinada) (Precio de escasez (VOLL))
    # ----------------------------------------------------------
    # 1-ter  (sin pendiente artificial; constante mantenida por compatibilidad)

    # ----------------------------------------------------------
    # 2. Construye automáticamente los coeficientes α_v→w
    #    a partir del CSV de tecnologías
    # ----------------------------------------------------------
    alpha = Dict{Tuple{Symbol,Symbol},Float64}()
    for tec in values(tecnologias)
        v_in  = Symbol(tec.vector_consumido)
        v_out = Symbol(tec.vector)
        @assert isfinite(tec.ratio_energia) "ratio_energia NaN en $(tec.nombre)"
        r = tec.ratio_energia
        @assert r > 0                  "ratio_energia ≤0 en $(tec.nombre)"
        if !isempty(tec.vector_consumido) && r > 0
            # Consumo interno coherente con el resto del código:
            #   consumo = producción / ratio_energia
            # ⇒ en el balance necesitamos α = 1/ratio_energia
            alpha[(v_in, v_out)] = 1.0 / r
        end
    end

    # ----------------------------------------------------------
    # 3. Modelo JuMP-Complementarity
    # ----------------------------------------------------------
    model = Model(HiGHS.Optimizer)
    set_silent(model)  # DESACTIVO LAS SALIDAS PORQUE YA FUNCIONA BIEN, ASÍ LIMPIO ESPACIO EN LA TERMINAL
    #set_optimizer_attribute(model, MOI.NumberOfThreads(), 8) # 8 es el máximo de hilos para HIGHS puesto así CREO QUE ES MÁS LENTO QUE SE CONFIGURE QUE EL SOLVER EN SÍ, TIEMPOS DE 0.5 s
    #set_optimizer_attribute(model, "parallel", "on") # activa paralelización
    # Comp. perfecta: solo variables de cantidad (coste social mínimo)


    demand_fixed = Dict{Symbol, Vector{Float64}}()     # DEMANDA INELÁSTICA
    for v in VECTORES
        # Para cada hora h en 0:23, obtenemos demanda o 0.0
        demand_fixed[v]  = [ get(demand_ext, (v, h), 0.0) for h in HORAS ]
        @assert all(isfinite, values(demand_fixed[v])) "Demanda del vector $v , hora $h no finita"

        #for h in HORAS println("Demanda del vector $v , hora $h = $(demand_fixed[v][h+1])") end
    end


    # 3-B  Clearing conjunto por vector-hora

    tec_out = Dict{String,Symbol}()
    tec_in  = Dict{String,Symbol}()
    for (t, tec) in tecnologias
        tec_out[t] = Symbol(tec.vector)
        tec_in[t]  = Symbol(tec.vector_consumido)
    end

    tecs_por_vector = Dict(v => String[] for v in VECTORES)
    for t in keys(tecnologias)
        push!(tecs_por_vector[tec_out[t]], t)
    end



    # 3-A  Variables de producción      q[a,t,h]
    @variable(model,
        0 <= q[a=AGENTES, t in keys(tecnologias), h=HORAS] <=
            cap[(a, Symbol(t), h)])


    # --------------------------------------------------
    # 3-A-2  Producción eléctrica total y síncrona por hora
    # --------------------------------------------------
    @expression(model, P_elec_total[h = HORAS],
        sum(q[a, t, h]
            for a in AGENTES,
            t in tecs_por_vector[:electricidad]))
            
    @expression(model, P_elec_sync[h = HORAS],
        sum(tecnologias[t].aporta_inercia * q[a, t, h]
            for a in AGENTES,
            t in tecs_por_vector[:electricidad]))


    # Busqueda de valores indebidos
    @assert all(isfinite, values(cap))      "cap contiene Inf/NaN"
    @assert all(isfinite, values(a_v))      "a_v contiene Inf/NaN"
    @assert all(v -> isfinite(v) && v>0, values(b_v))  "b_v contiene valores no finitos o ≤0"


    # 3-A-bis  Déficit físico (MWh no servidos)
    @variable(model, 0 <= d[v=VECTORES,h=HORAS] <= 1e8)



    # @constraint(model,
    #     cap_lim[a = AGENTES, t in keys(tecnologias), h = HORAS],
    #     q[a, t, h] <= cap[(a, Symbol(t), h)])
    # @constraint(model, cap_lim[a=AGENTES, t in keys(tecnologias), h=HORAS],
    #     q[a,t,h] <= cap[(a, Symbol(t), h)])



    # 3-A-1 Semillas iniciales 

    # Calcula la producción típica por vector-hora
    if day(fecha) == 1 && month(fecha) == 1
        prom_demanda = Dict(v => mean(demand_fixed[v]) for v in VECTORES)

        for a in AGENTES, t in keys(tecnologias), h in HORAS
            v = tec_out[t]
            # (a) cuota proporcional a la demanda fija de esa hora
            base = demand_fixed[v][h+1] / (length(AGENTES)*length(tecs_por_vector[v]))
            # (b) si la demanda es 0 → usa el 5 % de la media histórica
            fallback = 0.05 * prom_demanda[v]
            seed = max(base, fallback)            # nunca 0, pero “pequeño”
            # (c) respeta la capacidad real
            set_start_value(q[a,t,h], min(seed, cap[(a,Symbol(t),h)]))
        end

    else
        # ————————————————— WARM-START desde x0 (día anterior) —————————————————
        # EL ANTERIOR FUNCIONA BASTANTE BIEN, PERO AHORA QUE ES UN SISTEMA DE COMPETENCIA PERFECTA, SALVO EN CONTADAS OCASIONES, ES MEJOR UTILIZAR ESTE WARM-START
        #LAS ÚNICAS PARA LAS QUE PUEDE SER MEJOR EL OTRO ES EN CAMBIOS DE ESTACIONES O, SOBRE TOD, ANUALES
        if x0 !== nothing
            vars = all_variables(model)
            @assert length(vars) == length(x0) "La longitud de x0 ($(length(x0))) no coincide con nº de variables ($(length(vars)))"
            set_start_value.(vars, x0)   # JuMP/MOI: VariablePrimalStart
        else
            # Semilla por defecto (demanda promedio)
            set_start_value.(p, datos.demanda_promedio_p)
            set_start_value.(q, datos.demanda_promedio_q)
        end
        #——————————————————————————————————————————————————————————————————
    end

    # 3-B  Cantidad total neta Q
    @expression(model, Q[v = VECTORES, h = HORAS],
        # (+) Producción bruta del vector v
        sum(q[a,t,h] for a in AGENTES, t in tecs_por_vector[v])
        # (-) Consumo de v como insumo para otras tecnologías
        - sum(
            (alpha[(v, tec_out[t_consumidora])] * q[a, t_consumidora, h])
            for a in AGENTES for t_consumidora in keys(tecnologias)
            if tec_in[t_consumidora] == v && haskey(alpha, (v, tec_out[t_consumidora])); init=0.0)
    )
    

    # 3-C  Balance estricto oferta-demanda  (se activa d[v,h] si falta energía)
    @constraint(model, balance[v = VECTORES, h = HORAS],
        Q[v,h] + d[v,h] == demand_fixed[v][h+1])


    # --------------------------------------------------
    # 3-C-2  Restricción de inercia mínima (≥20% producción síncrona)
    # --------------------------------------------------
    @constraint(model, inercia_min[h = HORAS],
        P_elec_sync[h] ≥ MIN_SYNC_SHARE * P_elec_total[h])

    # 3-D  Función objetivo de coste social mínimo
    #       (coste marginal + multa por déficit)
    @objective(model, Min,
        sum(costo_marg[(a, Symbol(t))] * q[a,t,h]
            for a in AGENTES, t in keys(tecnologias), h in HORAS)
        +  sum(PENAL_DEF[v] * d[v,h] for v in VECTORES, h in HORAS))
    

    # ------------------------------------------------------------------
    # Curva de aprendizaje: aplica la reducción al coste marginal antes
    # de construir las condiciones KKT.
    # ------------------------------------------------------------------
    for tec in values(tecnologias)
        Tecnologias.actualizar_costo_marginal_aprendizaje!(tec)
    end

    # Sincroniza el diccionario de costes marginales que recibió la función
    # con los valores recién actualizados. El coste marginal es homogéneo
    # entre agentes (la curva depende de la tecnología), de modo que basta
    # con reasignarlo.
    for ag in AGENTES
        for (tec_nom, tec) in tecnologias
            costo_marg[(ag, Symbol(tec_nom))] = tec.costo_om_variable
        end
    end

    @info "Costes marginales actualizados (€/MWh): Min=$(minimum(values(costo_marg))), Max=$(maximum(values(costo_marg)))"

    ##################################################################
    # 3-E.  Resolución del problema social (HiGHS)
    ##################################################################
    optimize!(model)

    # Objeto devuelto por solveMCP (lo conservamos por si quieres inspeccionarlo)
    status = JuMP.termination_status(model)
    # Acceso a las cantidades óptimas:
    #    valor = value(q[a,t,h])
    # Si quieres guardarlas en un vector plano:
    #    q_opt = JuMP.value.(all_variables(model))

    #=status = solveMCP(model;           # ← NLsolve
    solver = :NLsolve,             # backend open-source
    show_trace   = Val(true),
    trace_level  = TraceAll(),
    autodiff     = AutoPolyesterForwardDiff(),
    nlsolve_options = (
        autodiff = :forward,
        convergence_tolerance = 10,  # mantiene tolerancia
        method        = :anderson,
        reformulation = :FischerBurmeister,       # NO (Fischer-Burmeister)
        xtol           = 10.05,      # ≈√(10²+1²): 10 MW + 1 €/MW
        ftol           = 1e-3,       # tolerancia para el residuo
        callback       = cb,
        )
    ) =#
    
    # ----------------------------------------------------------
    # 4. Precios resultantes
    # ----------------------------------------------------------
    if status == MOI.OPTIMAL
        P_final          = Dict{Tuple{Symbol,Int},Float64}()
        deficits_by_vect = Dict{String,Float64}()
        for v in VECTORES, h in HORAS
            # oferta total
            oferta = sum((value(q[a,t,h]) for a in AGENTES, t in tecs_por_vector[v]); init=0.0)
            # autoconsumo
            autoc  = sum(alpha[(tec_in[t],v)] * value(q[a,t,h])
                         for a in AGENTES, t in keys(tecnologias)
                         if haskey(alpha,(tec_in[t],v)); init=0.0)
            Q_net = oferta - autoc  # Prod neta
            # Precio marginal = dual de la restricción de balance
            P_final[(v,h)] = dual(balance[v,h]) * 1.1 # 10% margenes

            # déficit físico (MWh) si la oferta neta < demanda fija
            deficit = value(d[v,h])
            if deficit > 1e-3
                deficits_by_vect[string(v)] =
                    get(deficits_by_vect, string(v), 0.0) + deficit
            end
        end
        # ---- reconstruimos q (solo variables existentes en idx2info) ----
        q_dict = Dict{Tuple{Int,String,Int},Float64}()
        for a in AGENTES, t in keys(tecnologias), h in HORAS
            q_dict[(a, String(t), h)] = value(q[a,t,h])
        end
        # ---- construir despacho (agente, tecnología, vector, hora) ----
        dispatch_dict = Dict{Tuple{Int,String,String,Int},Float64}()
        for ((ag, tec, h), val) in q_dict
            vec_str = String(tec_out[tec])      # vector que produce esa tecnología
            dispatch_dict[(ag, tec, vec_str, h)] = val
        end

        return (
            q                 = q_dict,
            despacho          = dispatch_dict,
            precios_marginales = P_final,
            deficits          = deficits_by_vect,
            x0_next           = JuMP.value.(all_variables(model))
        )
    else
        error("El solver MCP no convergió")
    end
end

end # module