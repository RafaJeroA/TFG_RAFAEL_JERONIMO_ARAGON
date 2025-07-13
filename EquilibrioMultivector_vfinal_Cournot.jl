module EquilibrioMultivector
# Asegúrate de que NLsolve esté instalado:
# ] add NLsolve
using JuMP, Complementarity, NLsolve
#include("Datos.jl")
using ..Datos, ..Tecnologias
using Logging
using SparseArrays
using MathOptInterface
using SparseDiffTools


using Base.Threads

global_logger(ConsoleLogger(stderr, Logging.Debug))
const PENALTY_SLOPE_MULT = 1.0         # pendiente artificial desactivada
const FINE_DEFICIT  = 100.0         # €/MWh – multa base por déficit

const DEFICIT_TOL   = 10             # MWh  – tolerancia de déficit

const MOI = MathOptInterface

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
    model = MCPModel()

    # El precio deja de ser variable independiente: solo q permanece


    demand_fixed = Dict{Symbol, Vector{Float64}}()     # DEMANDA INELÁSTICA
    for v in VECTORES
        # Para cada hora h en 0:23, obtenemos demanda o 0.0
        demand_fixed[v]  = [ get(demand_ext, (v, h), 0.0) for h in HORAS ]
        @assert all(isfinite, values(demand_fixed[v])) "Demanda del vector $v , hora $h no finita"

        for h in HORAS println("Demanda del vector $v , hora $h = $(demand_fixed[v][h+1])") end
    end

    # 3-A  Restricción de capacidad + definición variable (cap)
    # q = prod[a,t,h]
    @variable(model,
    0 <= q[a=AGENTES, t in keys(tecnologias), h=HORAS]
      <= cap[(a, Symbol(t), h)])

    # Busqueda de valores indebidos
    @assert all(isfinite, values(cap))      "cap contiene Inf/NaN"
    @assert all(isfinite, values(a_v))      "a_v contiene Inf/NaN"
    @assert all(v -> isfinite(v) && v>0, values(b_v))  "b_v contiene valores no finitos o ≤0"

    # 3-A-1 Semillas iniciales 
    for a in AGENTES, t in keys(tecnologias), h in HORAS
        set_start_value(q[a,t,h], 0.1*cap[(a,Symbol(t),h)])
    end

    # 3-A-bis  Déficit físico (MWh no servidos)
    #@variable(model, d[v = VECTORES, h = HORAS] >= 0)
    # Techo = 150 % de la demanda fija o 20 000 MW si la demanda es 0
    @variable(model,
        0 <= d[v = VECTORES, h = HORAS] <=
            max(1.5 * demand_fixed[v][h+1], 20_000.0))


    # @constraint(model,
    #     cap_lim[a = AGENTES, t in keys(tecnologias), h = HORAS],
    #     q[a, t, h] <= cap[(a, Symbol(t), h)])
    # @constraint(model, cap_lim[a=AGENTES, t in keys(tecnologias), h=HORAS],
    #     q[a,t,h] <= cap[(a, Symbol(t), h)])

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

    
    # 3-B  Cantidad total Q y precio P(Q)
    @NLexpression(model, Q[v = VECTORES, h = HORAS],
        # (+) Producción bruta del vector v
        sum(q[a,t,h] for a in AGENTES, t in tecs_por_vector[v]; init=0.0)
        # (-) Consumo de v como insumo para otras tecnologías
        - sum(
            (alpha[(v, tec_out[t_consumidora])] * q[a, t_consumidora, h])
            for a in AGENTES for t_consumidora in keys(tecnologias)
            if tec_in[t_consumidora] == v && haskey(alpha, (v, tec_out[t_consumidora])); init=0.0)
    )
    

    # Saldo oferta-demanda
    @NLexpression(model, Δ[v = VECTORES, h = HORAS],
        Q[v,h] + d[v,h] - demand_fixed[v][h+1])

    # Precio inverso por tramos: lineal normal + lineal penalizado
    @NLexpression(model, P_inv[v = VECTORES, h = HORAS],
        a_v[v] - b_v[v]*Δ[v,h])
    

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

    # 3-C  Condiciones de optimalidad (KKT) de Cournot
    @mapping(model,
        F_q[a = AGENTES, t in keys(tecnologias), h = HORAS],
        # usamos el dict tec_out preparado más arriba para mapear t→vector
        -P_inv[tec_out[t], h]            # (–) precio: pasamos a minimización
        + b_v[tec_out[t]] * q[a,t,h]     # (+) efecto de mi producción en el precio
        + costo_marg[(a, Symbol(t))]     # (+) coste marginal (O&M + CO₂)
        + 1.0e-6 * q[a, t, h]            # (+) regularización numérica
    )
    @complementarity(model, F_q, q)


    # ─── Complementariedad para el déficit físico d ──────────────────────
    λ_pen = Dict((v,h) => FINE_DEFICIT for v in VECTORES, h in HORAS)

    @mapping(model,
        F_d[v = VECTORES, h = HORAS],
        λ_pen[(v,h)] - P_inv[v,h]        # ≥0 ⟂ d[v,h] ≥0
    )
    @complementarity(model,
        F_d,
        d
    )

    println("n variables  = ", length(all_variables(model)))
    println("restricciones no lineales: ", JuMP.num_nonlinear_constraints(model))

    cb = (state,) -> begin
        iter = state.iteration
        x    = state.zero_estimate      # vector actual de variables q y P empacado

        if any(!isfinite, x)
            @error "[NLsolve] NaN/Inf detectado en iter=$(iter); abortando"
            return true                       # se aborta la iteración
        end

        res  = norm(state.residual)     # norma del residuo complementario
        println("[NLsolve] iter=$iter  ‖residual‖=$(round(res, sigdigits=6))")
        return false                    # devuelve true si quieres abortar
    end


    ##################################################################
    # 3-D.  Resolución directa del MCP con NLsolve (robusta)
    ##################################################################

    for (k,v) in costo_marg
        @assert isfinite(v) "$k = $v no es finito"
    end
    #@assert all(isfinite, M)  "M contiene valores no finitos"
    #@assert all(isfinite, q_vec) "q_vec contiene valores no finitos"

    status = solveMCP(model;
        solver        = :NLsolve,
        autodiff = AutoSparseForwardDiff(),
        reformulation = :fischer,          # ← clave de estabilidad
        show_trace    = Val(true),
        nlsolve_options = (
            method   = :anderson,      # ← método más robusto
            ftol     = 1e-3,
            xtol     = 1e-3,
            callback = cb
        )
    )

    # Objeto devuelto por solveMCP (lo conservamos por si quieres inspeccionarlo)
    result = status

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
        convergence_tolerance = 10,  # mantiene tu tolerancia
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
    if JuMP.termination_status(model) == MOI.SUCCESS
        P_final          = Dict{Tuple{Symbol,Int},Float64}()
        deficits_by_vect = Dict{String,Float64}()
        for v in VECTORES, h in HORAS
            # oferta total
            oferta = sum((value(q[a,t,h]) for a in AGENTES, t in tecs_por_vector[v]); init=0.0)
            # autoconsumo
            autoc  = sum(alpha[(tec_in[t],v)] * value(q[a,t,h])
                         for a in AGENTES, t in keys(tecnologias)
                         if haskey(alpha,(tec_in[t],v)); init=0.0)
            Q_net = oferta - autoc
            Δ_vh  = Q_net - demand_fixed[v][h+1]
            P_final[(v,h)] = a_v[v] - b_v[v]*Δ_vh

            # déficit físico (MWh) si la oferta neta < demanda fija
            deficit = max(0.0, demand_fixed[v][h+1] - Q_net)
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