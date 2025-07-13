# =========================================
# OptimizacionSistema.jl - MILP multianual con JuMP
# =========================================

module OptimizacionSistema

# Dependencias externas
using JuMP, DataFrames, Statistics

using MathOptInterface
const MOI = MathOptInterface
#import MathOptInterface.Utilities as MOIU  
using  MathOptInterface:                              # los tipos viven en MOI
    ScalarAffineTerm, ScalarAffineFunction

import MultiObjectiveAlgorithms
const MOA = MultiObjectiveAlgorithms

# Número de agentes para repartir el coste
#const N_AGENTES = length(keys(agentes))

using MathOptInterface.Utilities: operate   


const CASH_MIN_OPER  = 30.0e6   # €
const RATIO_DEUDA_MAX = 3.5

function add_objectives_vector!(model::JuMP.Model,
        exprs::Vector{JuMP.AffExpr};
        sense = MOI.MAX_SENSE)
    comps = MOI.ScalarAffineFunction{Float64}[]
    for f in exprs
        # --- API nueva -----------------------------------------------------
        terms = [MOI.ScalarAffineTerm(a, JuMP.index(v))
            for (a, v) in JuMP.linear_terms(f)]
        # -------------------------------------------------------------------
        push!(comps, MOI.ScalarAffineFunction(terms, JuMP.constant(f)))
    end
    f_vec  = operate(vcat, Float64, comps...)
    backend = JuMP.backend(model)
    MOI.set(backend, MOI.ObjectiveSense(), sense)
    MOI.set(backend, MOI.ObjectiveFunction{typeof(f_vec)}(), f_vec)
end



#=
# --------------------------------------------------------------------
# Shim multi-objetivo para versiones antiguas de MOI (<1.14)
# --------------------------------------------------------------------
if !isdefined(MOI, :NumberOfObjectives)
    struct NumberOfObjectives <: MOI.AbstractModelAttribute end
    const _obj_counter = IdDict{Any,Int}()   # lleva la cuenta por modelo
    function _next_idx(model)
        _obj_counter[model] = get(_obj_counter, model, 0) + 1
    end
end

function add_objective(model::MOI.ModelLike,
                       sense::MOI.OptimizationSense,
                       f; weight::Real = 1.0,
                          priority::Integer = 1)
    idx = !isdefined(MOI, :NumberOfObjectives) ? _next_idx(model) :
          MOI.get(model, MOI.NumberOfObjectives()) + 1

    local f_moi
    if f isa JuMP.AffExpr                             # ----- lineal -----
        coeffs = [ScalarAffineTerm(a, JuMP.index(v))
                  for (a,v) in JuMP.linear_terms(f)]  # ⬅️ API moderna
        f_moi = ScalarAffineFunction(coeffs, JuMP.constant(f))  # constante
    elseif f isa JuMP.QuadExpr                    # Cuadrática (opcional)
        qterms = [ScalarQuadraticTerm(c, JuMP.index(v1), JuMP.index(v2))
                  for (c,v1,v2) in JuMP.quad_terms(f)]          # iterator
        aterms = [ScalarAffineTerm(a, JuMP.index(v))
                  for (a,v) in JuMP.linear_terms(f)]            # comp. lineal
        f_moi = ScalarQuadraticFunction(qterms, aterms, JuMP.constant(f))
    else                                          # Ya es MOI.AbstractFunction
        f_moi = f
    end
    MOI.set(model, MOI.ObjectiveFunction{typeof(f_moi)}(idx), f_moi)
    MOI.set(model, MOI.ObjectiveSense(idx),      sense)
    MOI.set(model, MOI.ObjectiveWeight(idx),     weight)
    MOI.set(model, MOI.ObjectivePriority(idx),   priority)
    return idx
end

# Versión que acepta `JuMP.Model`
function add_objective(model::JuMP.Model, args...; kw...)
    add_objective(JuMP.backend(model), args...; kw...)
end
# --------------------------------------------------------------------
=#

using CSV # Para exportar debug
using Base.Threads # Para nthreads() y @threads
using ..Agentes  # para cargar estructuras y parámetros de agentes

# Dependencias internas del proyecto (usando .)
using ..Tecnologias # Acceder a la estructura Tecnologia
using ..Escenarios # Acceder a la estructura Escenario
using ..InversionesPendientes # Para PENDING_INV
using Dates
using Main.Utils: normalizar_nombre_tecnologia
using ..FinanzasParametros: EQ_SHARE, INT_RATE, PLAZO_PRESTAMO


# Solver (asegúrate que HiGHS está instalado en el environment)
import HiGHS
using HiGHS
const MAX_SIMPLEX_CONCURRENCY = 8  # límite interno de HiGHS


export ejecutar_optimizacion_agente, ejecutar_optimizacion_multiagente, construir_submodelo_agente!
export acumular_gen_verde!


# ─── Capacidad existente a nivel sistema ────────────────────────────────
# Diccionario global: tecnología → Dict{Periodo,Float64}
const Periodo = Tuple{Int,Int}
const cap_existente = Dict{String,Dict{Periodo,Float64}}()


const MAX_RATIO_DEUDA = 3.5      # ≤ 3.5× patrimonio
const DSCR_MAX        = 0.40     # ≤ 40 % ingresos
#const CAP_MIN_MW      = 0     # ≥ 10 MW

# ───────── Límites de inversión (replica heurística del sub-MILP) ─────────
const LIM_MIN_MW      = 100.0      # techo absoluto mensual (MW)
const LIM_FRAC_POT    = 0.03       # 3 % de la potencia propia existente
const LIM_FRAC_CASH   = 0.20       # máx 20 % del cash disponible
const PENAL_EXCESO    = 1.0e3     # € por MW ó € de cash fuera de límite

# ─────────  Variable de demanda no servida ─────────
const COSTE_DNS = 1.0e5        # €/MWh; ajusta si lo deseas
const COSTE_SOBREPRODUCCION = 1.0e3        # €/MWh; ajusta si lo deseas


export get_cap_existente, get_cap_existente_agente
"""
get_cap_existente(tecnologia::String, periodo::Periodo) → Float64
Devuelve la capacidad total instalada de `tecnologia` en el mes `periodo`.
"""
get_cap_existente(tecnologia::String, periodo::Periodo) =
get(
get(cap_existente, tecnologia, Dict{Periodo,Float64}()),
periodo,
0.0
)

"""
get_cap_existente(tecnologia::String, anio::Int, mes::Int) → Float64
Sobrecarga que acepta año y mes por separado.
"""
get_cap_existente(tecnologia::String, anio::Int, mes::Int) =
get_cap_existente(tecnologia, (anio, mes))

# ─── Capacidad existente por agente ─────────────────────────────────────
"""
get_cap_existente_agente(ag::Int, tecnologia::String, periodo::Periodo, agentes::Dict{Int,Agentes.EnergyAgent}) → Float64
Devuelve la capacidad de `tecnologia` controlada por el agente `ag` en `periodo`.
"""
function get_cap_existente_agente(id::Int, tec::String, periodo::Periodo,
    agentes::Dict{Int,Agentes.EnergyAgent})
    cap = get(agentes[id].capacidades, tec, 0.0)
    return cap isa Number ? cap : get(cap, periodo, 0.0)
end

function ejecutar_optimizacion_agente(
    agente_obj::Main.Agentes.EnergyAgent,
    datos_agente_especificos::NamedTuple,
    datos_globales_bundle,
    tecnologias_del_agente::Dict{String,Main.Tecnologias.Tecnologia},
    escenario_actual::Main.Escenarios.Escenario,
    precios_mercado_esperados::Dict{Tuple{Int,String,String,Int},Float64};
    tiempo_max::Int = 600)

    model = Model(() -> MOA.Optimizer(HiGHS.Optimizer))
    # ───────── Configuración de HiGHS ─────────
    set_optimizer_attribute(model, "time_limit", float(tiempo_max))
    set_optimizer_attribute(model, "threads", 1)       # o Threads.nthreads()
    #set_optimizer_attribute(model, "parallel", "on")   # habilita paralelismo interno
    set_optimizer_attribute(model, "solver", "simplex")           
    #set_optimizer_attribute(model, "simplex_strategy", "dual")
    set_optimizer_attribute(model, "run_crossover", "off")   # sin crossover
    set_optimizer_attribute(model, "threads", 1)              # UN solo hilo
    set_optimizer_attribute(model, "parallel", "off")        # desactiva paralelismo
    set_silent(model)                                        # sin logs del solver

    info = construir_submodelo_agente!(model, agente_obj, datos_agente_especificos,
        datos_globales_bundle, tecnologias_del_agente, escenario_actual,
        precios_mercado_esperados, "ag$(agente_obj.id)_", agentes)

    # Añadimos único objetivo y guardamos su índice (para poder leerlo luego)
    add_objectives_vector!(model, [info.beneficio])

    optimize!(model)

    status_final = termination_status(model)
    obj_final    = has_values(model) ? objective_value(model) : NaN

    cap_df, prod_df, alm_df = construir_resultados_agente(info.id_agente,
        info.tecs_op, info.tecs_alm, info.PERIODOS, info.T_global,
        info.cap_total, info.cap_nueva,
        info.produccion, info.carga, info.descarga, info.nivel)

    return (
        agente_id     = info.id_agente,
        escenario     = escenario_actual.nombre,
        opt_status    = status_final,
        opt_objective = obj_final, # MILLONES €
        opt_resultados = Dict(:capacidad=>cap_df,:produccion=>prod_df,:almacenamiento=>alm_df),
        model = model,
    )
end

function ejecutar_optimizacion_multiagente(
    agentes::Dict{Int,Main.Agentes.EnergyAgent},
    datos_por_agente::Dict{Int,T} where {T<:NamedTuple},
    datos_globales_bundle,
    tecnologias_por_agente::Dict{Int,Dict{String,Main.Tecnologias.Tecnologia}},
    escenario_actual::Main.Escenarios.Escenario,
    precios_mercado_esperados::Dict{Tuple{Int,String,String,Int},Float64};
    tiempo_max::Int = 3600)


    model = Model(() -> MOA.Optimizer(HiGHS.Optimizer))
    # ───────── Configuración de HiGHS ─────────
    set_optimizer_attribute(model, "time_limit", float(tiempo_max))
    set_optimizer_attribute(model, "threads", 1)       # o Threads.nthreads()
    #set_optimizer_attribute(model, "parallel", "on")   # habilita paralelismo interno
    set_optimizer_attribute(model, "solver", "simplex")           
    #set_optimizer_attribute(model, "simplex_strategy", "dual")

    set_optimizer_attribute(model, "run_crossover", "off")   # desactivamos crossover
    set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "parallel", "off")        # desactiva paralelismo interno

    sub       = Dict{Int,NamedTuple}()
    #objetivos = JuMP.AffExpr[]
    for id in sort(collect(keys(agentes)))
        info = construir_submodelo_agente!(model, agentes[id], datos_por_agente[id],
            datos_globales_bundle, tecnologias_por_agente[id], escenario_actual,
            precios_mercado_esperados, "ag$(id)_", agentes)
        #push!(objetivos, info.beneficio)
        sub[id] = info
    end

    demanda_df = datos_globales_bundle.demanda
    demanda_map = Dict{Tuple{Int,String,String,Int},Float64}()
    for r in eachrow(demanda_df)
        demanda_map[(Int(r.anio), String(r.vector), String(r.estacion), Int(r.hora))] = Float64(r.valor)
    end

    PERIODOS = first(values(sub)).PERIODOS
    T_global = first(values(sub)).T_global
    vectores = unique(String.(demanda_df.vector))


    # ─────────  Variable de demanda no servida ─────────
    @variable(model,
        dns[p in PERIODOS,
            vec in vectores,
            (e,h) in T_global] >= 0,
            base_name = "DNS")

    @variable(model,
        sobreproduccion[p in PERIODOS,
            vec in vectores,
            (e,h) in T_global] >= 0,
            base_name = "Sobreproduccion")


    # 1) Cuota de capacidad de cada agente (excluye import_*)
    # @expression(model,
    #     cap_total_sin_importes[ag in keys(sub), (a,m) in PERIODOS],
    #         (sum(sub[ag].cap_total[tec, (a,m)]
    #             for tec in sub[ag].tecs_op        # todas las tecs del agente
    #             if !startswith(tec,"import_")); init = 0.0)    # excluye import_*
    # )

    # 2. Proporción sobre el total del sistema
    #@expression(model,
    #    propor_cap[ag in keys(sub), (a,m) in PERIODOS],
    #        cap_total_sin_importes[ag,(a,m)] /
    #        sum(cap_total_sin_importes[id,(a,m)] for id in keys(sub))
    #)

    # 3) Penalización individual de la DNS
    @expression(model,
        penalizacion_dns_ag[ag in keys(sub)],
            sum(
                COSTE_DNS/1e6 * dns[p, vec, (e,h)] / 24
                for p in PERIODOS, vec in vectores, (e,h) in T_global
            )
            #sum(propor_cap[ag,(a,m)] * dns[(a,m), vec, (e,h)]
            #    for (a,m) in PERIODOS, vec in vectores, (e,h) in T_global )
    )

    @expression(model,
    penalizacion_sobreproduccion_ag[ag in keys(sub)],
        sum(
            COSTE_SOBREPRODUCCION/1e6 * sobreproduccion[p, vec, (e,h)] / 20
            for p in PERIODOS, vec in vectores, (e,h) in T_global
        )
    )

    # 4) Nuevo vector de objetivos (beneficio ya con la penalización)
    objetivos_final = [
        sub[id].beneficio - penalizacion_dns_ag[id] - penalizacion_sobreproduccion_ag[id]
        for id in sort(collect(keys(sub)))
    ]

    add_objectives_vector!(model, objetivos_final)


    # ─────────  Penalización de demanda no servida ─────────
    #@expression(model, penalizacion_dns,
    #    -COSTE_DNS *
    #    sum(dns[(a,m), vec, (e,h)]
    #        for (a,m) in PERIODOS, vec in vectores, (e,h) in T_global))

    #@expression(model, beneficio_total, sum(objetivos))
    #@objective(model, Max, beneficio_total)
    #push!(objetivos, penalizacion_dns)
    
    # Finalmente, crea el vector multi-objetivo
    #add_objectives_vector!(model, objetivos)


    for p in PERIODOS, vec in vectores, (e,h) in T_global
        (a,m) = p
        demanda = get(demanda_map, (a, vec, e, h), 0.0)
        if demanda > 0
            @constraint(model,
                        # ① Generación directa de todas las tecnologías cuyo vector coincide
                        sum(info.produccion[t,p,(e,h)]
                            for (id, info) in sub, t in info.tecs_op
                            if tecnologias_por_agente[id][t].vector == vec)
            
                        # ② Flujo neto de almacenamiento (descarga − carga)
                      + sum(info.descarga[t,p,(e,h)] - info.carga[t,p,(e,h)]
                            for (id, info) in sub, t in info.tecs_alm
                            if tecnologias_por_agente[id][t].vector == vec)
            
                        # ③ Demanda no servida (variable recién creada)
                      + dns[p, vec, (e,h)]
            
                    == demanda + sobreproduccion[p, vec, (e,h)])
        end
    end

    optimize!(model)

    let
      ids = sort(collect(keys(sub)))
      z   = JuMP.objective_value(model)        # ahora z[i] corresponde a ids[i]
      for (pos, id) in enumerate(ids)
        info = sub[id]
        cap_df, prod_df, alm_df = construir_resultados_agente(
          id, info.tecs_op, info.tecs_alm,
          info.PERIODOS, info.T_global,
          info.cap_total, info.cap_nueva,
          info.produccion, info.carga, info.descarga, info.nivel
        )
        sub[id] = merge(info, (
          opt_status    = termination_status(model),
          opt_objective = z[pos],            # ✔️ usar pos, no el ID
          opt_resultados = Dict(
            :capacidad     => cap_df,
            :produccion    => prod_df,
            :almacenamiento=> alm_df
          ),
        ))
      end
    end

            dns_df = DataFrame(anio = Int[], mes = Int[], vector = String[],
            estacion = String[], hora = Int[], dns_mw = Float64[])
        for p in PERIODOS, v in vectores, (e,h) in T_global
        (a,m) = p
        push!(dns_df, (a, m, v, e, h, value(dns[p, v, (e,h)])))
        end
        CSV.write("dns_resultados.csv", dns_df)

        sobreproduccion_df = DataFrame(anio = Int[], mes = Int[], vector = String[],
            estacion = String[], hora = Int[], sobreproduccion_mw = Float64[])
        for p in PERIODOS, v in vectores, (e,h) in T_global
        (a,m) = p
        push!(sobreproduccion_df, (a, m, v, e, h, value(sobreproduccion[p,v,(e,h)])))
        end
        CSV.write("sobreproduccion_resultados.csv", sobreproduccion_df)


    return (model=model, sub=sub)
end


"""
get_cap_existente_agente(ag::Int, tec::String, anio::Int, mes::Int, agentes::Dict{Int,Agentes.EnergyAgent}) → Float64
Sobrecarga que acepta año y mes por separado.
"""
get_cap_existente_agente(ag::Int, tec::String, anio::Int, mes::Int, agentes::Dict{Int,Agentes.EnergyAgent}) =
get_cap_existente_agente(ag, tec, (anio, mes), agentes)



"""
construir_submodelo_agente!(
    model::JuMP.Model,
    agente_obj::Main.Agentes.EnergyAgent,
    datos_agente_especificos::NamedTuple,
    datos_globales_bundle,
    tecnologias_del_agente::Dict{String,Main.Tecnologias.Tecnologia},
    escenario_actual::Main.Escenarios.Escenario,
    precios_mercado_esperados::Dict{Tuple{Int,String,String,Int},Float64},
    prefijo::String,
)
Construye variables y restricciones del MILP para un único agente dentro de
`model`. Los nombres de las variables se prefijan para garantizar unicidad.
"""
function construir_submodelo_agente!(
    model::JuMP.Model,
    agente_obj::Main.Agentes.EnergyAgent,
    datos_agente_especificos::NamedTuple,
    datos_globales_bundle,
    tecnologias_del_agente::Dict{String, Main.Tecnologias.Tecnologia},
    escenario_actual::Main.Escenarios.Escenario,
    precios_mercado_esperados::Dict{Tuple{Int, String, String, Int}, Float64},
    prefijo::String,
    agentes::Dict{Int,Main.Agentes.EnergyAgent}
)
    nombre_agente = agente_obj.name
    id_agente = agente_obj.id
    cap_existente_agente_actual = datos_agente_especificos.cap_existente_agente

    # ─── Función local para leer capacidad histórica del propio agente ──────
    get_cap_existente_agente(id, tec, anio, mes) =
        get(
            get(cap_existente_agente_actual, tec, Dict{Tuple{Int,Int},Float64}()),
            (anio, mes),
            0.0
        )

    println("🧩 Construyendo MILP para maximizar beneficio del AGENTE: $id_agente - $nombre_agente en escenario: $(escenario_actual.nombre)...")

    demanda_total_sistema_df = datos_globales_bundle.demanda

    # Periodos de simulación: tuplas (año, mes)
    PERIODOS = sort(
        unique(
          [(Int(r.anio), Int(r.mes)) for r in eachrow(demanda_total_sistema_df)]
        )
      )

    perfiles_gen_df_global = datos_globales_bundle.perfiles_generacion
    politicas_df_global = datos_globales_bundle.politicas
    tasa_descuento = datos_globales_bundle.tasa_descuento

    #ANIOS = sort(unique(Int.(demanda_total_sistema_df.anio)))
    TECNOLOGIAS_OPERABLES_AGENTE = collect(keys(tecnologias_del_agente))

    # ─── BLOQUEAR NUEVAS INVERSIONES EN TECNOLOGÍAS CON RETIRO PROGRAMADO ───
    if hasproperty(escenario_actual, :retiro_programado)
        tecs_retiro = keys(escenario_actual.retiro_programado)
        if !isempty(tecs_retiro)
            println("    🔒 Bloqueando inversiones en: $(join(tecs_retiro, ", ")) para $(agente_obj.name)")
            filter!(tec -> !(tec in tecs_retiro), TECNOLOGIAS_OPERABLES_AGENTE)
        end
    end
    # ─── FIN BLOQUEO ─────────────────────────────────────────────────────────

    # ─── Diccionarios con límites mensuales (MW y €) ──────────────────────────
    lim_cap  = Dict{Tuple{String,Periodo},Float64}()   # (tec,(a,m)) ⇒ MW
    lim_cash = Dict{Periodo,Float64}()                 # (a,m)      ⇒ €

    for p in PERIODOS
        (a,m) = p
        lim_cash[p] = agente_obj.cash * LIM_FRAC_CASH
        for t in TECNOLOGIAS_OPERABLES_AGENTE
            # capacidad instalada propia al final del mes anterior
            cap_prev = (m == 1) ?
                get_cap_existente_agente(id_agente, t, a-1, 12) :
                get_cap_existente_agente(id_agente, t, a, m-1)

            lim_cap[(t,p)] = max(LIM_FRAC_POT * cap_prev, LIM_MIN_MW)
        end
    end

    if isempty(TECNOLOGIAS_OPERABLES_AGENTE)
        @warn "Agente $nombre_agente no tiene tecnologías operables. Saltando optimización."
        return (agente_id=id_agente, escenario=escenario_actual.nombre, opt_status=:NO_TECNOLOGIAS, opt_objective=NaN, opt_resultados=Dict{Symbol,DataFrame}(), model=nothing)
    end

    # Determinar ESTACIONES a partir de la demanda global o pasarlo en el bundle si es más robusto
    # Esta línea asume que Datos.ESTACIONES está disponible y es correcto
    ESTACIONES_DISPONIBLES = unique(String.(demanda_total_sistema_df.estacion))
    T_global = unique([(String(r.estacion), Int(r.hora)) for r in eachrow(demanda_total_sistema_df)])

    demanda_dict_global_sys = Dict{Tuple{Int, String, String, Int}, Float64}()
    for r_dem in eachrow(demanda_total_sistema_df)
        # Clave: (anio, vector, estacion, hora)
        key_dem = (Int(r_dem.anio), String(r_dem.vector), String(r_dem.estacion), Int(r_dem.hora))
        demanda_dict_global_sys[key_dem] = Float64(r_dem.valor)
    end
    get_demanda_sistema(a, vec, e, h, def=0.0) = get(demanda_dict_global_sys, (a, vec, e, h), def)

    


    # dentro de ejecutar_optimizacion_agente, justo después de:
    perfiles_gen_df_global = datos_globales_bundle.perfiles_generacion

    perfiles_gen_global_map = Dict{Tuple{String, String}, Dict{Int, Float64}}()

    # 1) definimos las columnas requeridas
    req_perf_cols  = [:tecnologia, :estacion, :hora, :factor]
    # 2) convertimos los names() (que son Strings) a Symbols
    perf_cols_sym  = Symbol.(names(perfiles_gen_df_global))

    # 3) comprobación
    if !all(col -> col in perf_cols_sym, req_perf_cols)
        @warn "Faltan columnas requeridas en perfiles_gen_df_global: " *
            "$(setdiff(req_perf_cols, perf_cols_sym)) para el agente $nombre_agente."
    else
        # 4) agrupar y poblar el Dict
        for group_perf in groupby(perfiles_gen_df_global, [:tecnologia, :estacion])
            tecnologia = String(first(group_perf.tecnologia))
            est = String(first(group_perf.estacion))
            perfiles_gen_global_map[(tecnologia, est)] =
                Dict(Int(row.hora) => Float64(row.factor) for row in eachrow(group_perf))
        end
    end




    # ─── Precio de CO₂ por año ───────────────────────────────────────────────
    precio_co2_map = Dict{Int,Float64}()

    if :precio_co2 in propertynames(politicas_df_global) &&
    :anio      in propertynames(politicas_df_global)

        for r in eachrow(politicas_df_global)
            precio_co2_map[Int(r.anio)] = Float64(r.precio_co2)
        end
    else
        @warn "No se encontró la columna :precio_co2 en politicas.csv; se usa 0 €/t."
        for (a, m) in PERIODOS
            precio_co2_map[a] = 0.0
        end
    end
    # ─────────────────────────────────────────────────────────────────────────

    get_precio_co2_val(a::Int, def=0.0) = get(precio_co2_map, a, def)
    
  
    
    # ──────────────────────────────────────────────
    # Variables del agente *sin* registrar el nombre
    # (evita colisiones cuando se repite el bloque)
    # ──────────────────────────────────────────────
    cap_nueva_ag = @variable(model,
        [tecnologia in TECNOLOGIAS_OPERABLES_AGENTE,
         p      in PERIODOS],
        lower_bound = 0,
        base_name   = "$(prefijo)cap_nueva_")
    
    cap_total_ag = @variable(model,
        [tecnologia in TECNOLOGIAS_OPERABLES_AGENTE,
         p      in PERIODOS],
        lower_bound = 0,
        base_name   = "$(prefijo)cap_total_")
    
    produccion_ag = @variable(model,
        [tecnologia in TECNOLOGIAS_OPERABLES_AGENTE,
         p      in PERIODOS,
         (e,h)      in T_global],
        lower_bound = 0,
        base_name = "$(prefijo)produccion_")
    
    exceso_cap = @variable(model,
        [t        in TECNOLOGIAS_OPERABLES_AGENTE,
         p    in PERIODOS],
        lower_bound = 0,
        base_name = "$(prefijo)exceso_cap_")
    
    exceso_cash = @variable(model,
        [p in PERIODOS],
        lower_bound = 0,
        base_name = "$(prefijo)exceso_cash_")

    # ───── Restricción de *límite mensual* (4 % o 100 MW) ─────
    let p0 = first(PERIODOS)
        (a0,m0) = p0
        for tec in TECNOLOGIAS_OPERABLES_AGENTE
            @constraint(model, cap_nueva_ag[tec,p0] <= max(0.04 * get_cap_existente_agente(id_agente, tec, a0, m0), 100.0))
        end
    end
    # ───── Horizonte rodante: bloquear inversiones fuera del 1er periodo ─────
    primer_periodo = first(PERIODOS)
    for p in PERIODOS
        if p != primer_periodo
            for t in TECNOLOGIAS_OPERABLES_AGENTE
                @constraint(model, cap_nueva_ag[t,p] == 0)
            end
        end
    end
        


    # ─── Límites de inversión mensuales ───────────────────────────────────────
    @constraint(model,
        [t in TECNOLOGIAS_OPERABLES_AGENTE, p in PERIODOS],
        cap_nueva_ag[t,p] <= lim_cap[(t,p)] + exceso_cap[t,p])

    @constraint(model,
        [p in PERIODOS],
        sum(cap_nueva_ag[t,p] * tecnologias_del_agente[t].costo_inicial
            for t in TECNOLOGIAS_OPERABLES_AGENTE) <=
        lim_cash[p] + exceso_cash[p])


    # ─── Expresiones financieras ─────────────────────────────────────────
    deuda_nueva_ag = @expression(model,
        [p in PERIODOS],
        sum((1 - EQ_SHARE) *
            tecnologias_del_agente[tec].costo_inicial *
            cap_nueva_ag[tec, p] for tec in TECNOLOGIAS_OPERABLES_AGENTE))

    deuda_total_ag = @expression(model,
        [p in PERIODOS],
        agente_obj.debt + deuda_nueva_ag[p])
    
    cuota_nueva_ag = @expression(model,
        [p in PERIODOS],
        deuda_nueva_ag[p] * (INT_RATE*(1+INT_RATE)^PLAZO_PRESTAMO) / ((1+INT_RATE)^PLAZO_PRESTAMO - 1))

    cuota_existente_ag = @expression(model,
        [p in PERIODOS],
        agente_obj.debt * (INT_RATE*(1+INT_RATE)^PLAZO_PRESTAMO) 
                          / ((1+INT_RATE)^PLAZO_PRESTAMO - 1))

    cuota_total_ag = @expression(model,
        [p in PERIODOS],
        cuota_nueva_ag[p] + cuota_existente_ag[p])

    cap_total_mw_ag = @expression(model,
        [p in PERIODOS],
        sum(cap_total_ag[tec, p] for tec in TECNOLOGIAS_OPERABLES_AGENTE))


    cap_total_sin_importes = @expression(model,
        [p in PERIODOS],
            (sum(cap_total_ag[tec, p]
                for tec in TECNOLOGIAS_OPERABLES_AGENTE
                if !startswith(tec,"import_")); init = 0.0))    # excluye import_*

    patrimonio_neto_ag = @expression(model, [p in PERIODOS],
        cap_total_sin_importes[p] * 7.8894e6 + (cap_total_mw_ag[p] - cap_total_sin_importes[p]) * 1000)   # 60 €/MWh × 24 h × 365 d × 15 a
    ##########################################################################

    TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE = [tecnologia for tecnologia in TECNOLOGIAS_OPERABLES_AGENTE if tecnologias_del_agente[tecnologia].tipo == "almacenamiento"]
    if !isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
        carga_ag = @variable(model,
            [tecnologia in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
             p in PERIODOS,
             (e,h)      in T_global],
            lower_bound = 0,
            base_name = "$(prefijo)carga_")
        
        descarga_ag = @variable(model,
            [tecnologia in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
             p in PERIODOS,
             (e,h)      in T_global],
            lower_bound = 0,
            base_name = "$(prefijo)descarga_")
        nivel_ag = @variable(model,
            [tecnologia in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
             p in PERIODOS,
             (e,h)      in T_global],
            lower_bound = 0,
            base_name = "$(prefijo)nivel_")
    end
    
    # Restricciones del agente
    for tec_ag_bal in TECNOLOGIAS_OPERABLES_AGENTE, p in PERIODOS
        (anio_bal, mes_bal) = p
        vida_util_tec_ag = max(1, tecnologias_del_agente[tec_ag_bal].vida_util)
        idx_bal = findfirst(==( (anio_bal, mes_bal) ), PERIODOS)
        if idx_bal == 1
            cap_inicial_val_ag = get_cap_existente_agente(id_agente, tec_ag_bal, anio_bal, mes_bal)
            @constraint(model, cap_total_ag[tec_ag_bal, (anio_bal, mes_bal)] == cap_inicial_val_ag + cap_nueva_ag[tec_ag_bal, (anio_bal, mes_bal)])
        else
            cap_prev_val_ag = cap_total_ag[tec_ag_bal, PERIODOS[idx_bal-1]]
            anio_retiro_val_ag = anio_bal - vida_util_tec_ag
            retiro_val_ag = (anio_retiro_val_ag >= first(PERIODOS)[1]) ? cap_nueva_ag[tec_ag_bal, (anio_retiro_val_ag, mes_bal)] : 0.0
            
            cap_inicial_original_ag = get_cap_existente_agente(id_agente, tec_ag_bal, anio_bal, mes_bal)
            retiro_inicial_especifico_ag = 0.0
            if anio_retiro_val_ag < first(PERIODOS)[1] && (first(PERIODOS)[1] + vida_util_tec_ag == anio_bal)
                 retiro_inicial_especifico_ag = cap_inicial_original_ag
            end
            @constraint(model, cap_total_ag[tec_ag_bal, (anio_bal, mes_bal)] == cap_prev_val_ag + cap_nueva_ag[tec_ag_bal, (anio_bal, mes_bal)] - retiro_val_ag - retiro_inicial_especifico_ag)
        end
    end

    if !isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
        @constraint(model,
            [tec in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
             p in PERIODOS,
             (e,h) in T_global],
            produccion_ag[tec,p,(e,h)] == 0)
    end
    


    for tec_prod_ag in TECNOLOGIAS_OPERABLES_AGENTE, (anio_prod, mes_prod) in PERIODOS, (e_prod,h_prod) in T_global
        if !(tec_prod_ag in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
            perfil_tec_ag = get(perfiles_gen_global_map, (tec_prod_ag, e_prod), nothing)
            factor_prod_ag = isnothing(perfil_tec_ag) ? 1.0 : get(perfil_tec_ag, h_prod, 0.0)
            @constraint(model, produccion_ag[tec_prod_ag, (anio_prod, mes_prod), (e_prod,h_prod)] <= cap_total_ag[tec_prod_ag, (anio_prod, mes_prod)] * factor_prod_ag)
        end
    end
    
    if !isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
        T_sorted_ag_alm = sort(T_global, by = x -> (findfirst(==(x[1]), ESTACIONES_DISPONIBLES), x[2]))
        for tec_alm_ag in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE
            eficiencia_ag_alm  = clamp(tecnologias_del_agente[tec_alm_ag].eficiencia, 0.01, 1.0) # Evitar división por cero si eficiencia es 0
            eff_sqrt_ag_alm    = sqrt(eficiencia_ag_alm)

            cap_ratio_ag_alm = 10.0        # 10 MWh por MW instalado

            cap_inicial_val_ag_alm = get_cap_existente_agente(id_agente, tec_alm_ag, first(PERIODOS)[1], first(PERIODOS)[2])

            for (anio, mes) in PERIODOS
                # Condición inicial de nivel al principio de cada año
                idx_period = findfirst(==( (anio, mes) ), PERIODOS)
                @constraint(model, nivel_ag[tec_alm_ag, (anio, mes), T_sorted_ag_alm[1]] ==
                    (if (anio,mes) == first(PERIODOS)
                            0.5 * cap_inicial_val_ag_alm * cap_ratio_ag_alm # Estado inicial para el primer año
                    else
                        nivel_ag[tec_alm_ag, PERIODOS[idx_period-1], T_sorted_ag_alm[end]] * (1 - tecnologias_del_agente[tec_alm_ag].tasa_autodescarga) # Nivel final del periodo anterior
                    end)
                    + carga_ag[tec_alm_ag, (anio, mes), T_sorted_ag_alm[1]] * eff_sqrt_ag_alm
                    - descarga_ag[tec_alm_ag, (anio, mes), T_sorted_ag_alm[1]] / eff_sqrt_ag_alm
                )

                for idx_alm_ag in 2:length(T_sorted_ag_alm) # Empezar desde el segundo paso temporal del año
                    (e_curr, h_curr) = T_sorted_ag_alm[idx_alm_ag]
                    (e_prev_alm, h_prev_alm) = T_sorted_ag_alm[idx_alm_ag-1]
                    nivel_prev_val_ag = nivel_ag[tec_alm_ag, (anio, mes), (e_prev_alm, h_prev_alm)]
                    
                    # ESTE ES ANTIGUO @constraint(model, nivel_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] == nivel_prev_val_ag + carga_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] * eff_sqrt_ag_alm - descarga_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] / eff_sqrt_ag_alm)
                    @constraint(model,
                        nivel_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] ==
                        nivel_prev_val_ag * (1 - tecnologias_del_agente[tec_alm_ag].tasa_autodescarga)
                      + carga_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] * eff_sqrt_ag_alm
                      - descarga_ag[tec_alm_ag, (anio, mes), (e_curr,h_curr)] / eff_sqrt_ag_alm
                    )
                end
                # Restricciones comunes para todos los pasos temporales del año
                for (e_alm_loop, h_alm_loop) in T_sorted_ag_alm
                    ctot_ag_alm = cap_total_ag[tec_alm_ag, (anio, mes)]
                    @constraint(model, nivel_ag[tec_alm_ag, (anio, mes), (e_alm_loop,h_alm_loop)] <= ctot_ag_alm * cap_ratio_ag_alm)
                    @constraint(model, carga_ag[tec_alm_ag, (anio, mes), (e_alm_loop,h_alm_loop)] <= ctot_ag_alm)
                    @constraint(model, descarga_ag[tec_alm_ag, (anio, mes), (e_alm_loop,h_alm_loop)] <= ctot_ag_alm)
                end
            end
        end
    end

    for vec_r in unique(t.vector for t in values(tecnologias_del_agente)), (anio_r, mes_r) in PERIODOS, (e_r,h_r) in T_global
        produccion_total_agente_vector = sum(produccion_ag[tecnologia,(anio_r,mes_r),(e_r,h_r)] for tecnologia in TECNOLOGIAS_OPERABLES_AGENTE if tecnologias_del_agente[tecnologia].vector == vec_r && !(tecnologia in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE); init=0.0)
        descarga_total_agente_vector = 0.0
        if !isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
            descarga_total_agente_vector = sum(descarga_ag[tecnologia,(anio_r,mes_r),(e_r,h_r)] for tecnologia in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE if tecnologias_del_agente[tecnologia].vector == vec_r; init=0.0)
        end




        # NUEVO: energía absorbida al cargar almacenamiento
        carga_total_agente_vector = isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) ? 0.0 :
            sum(carga_ag[tec, (anio_r, mes_r), (e_r,h_r)]
                for tec in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE
                if tecnologias_del_agente[tec].vector == vec_r ; init = 0.0)







        demanda_vector_sistema = get_demanda_sistema(anio_r, vec_r, e_r, h_r)
        # Para el punto 2, el agente no debería superar la demanda total del sistema.
        # Para el punto 3 (Cournot), esta 'demanda_vector_sistema' sería la demanda *residual* para el agente.
        @constraint(model,
                produccion_total_agente_vector + descarga_total_agente_vector
                    <= demanda_vector_sistema + carga_total_agente_vector)
    end
    
    coste_inv_ag = @expression(model, [p in PERIODOS], sum(cap_nueva_ag[t,p]*tecnologias_del_agente[t].costo_inicial for t in TECNOLOGIAS_OPERABLES_AGENTE; init=0.0))
    coste_omf_ag = @expression(model, [p in PERIODOS], sum(cap_total_ag[t,p]*tecnologias_del_agente[t].costo_om for t in TECNOLOGIAS_OPERABLES_AGENTE; init=0.0))
    
    costo_om_var_expr_terms = AffExpr()
    for (a_omv, m_omv) in PERIODOS
        for t_omv in TECNOLOGIAS_OPERABLES_AGENTE
            if !(t_omv in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
                for (e_omv,h_omv) in T_global
                    add_to_expression!(costo_om_var_expr_terms, produccion_ag[t_omv,(a_omv,m_omv),(e_omv,h_omv)] * tecnologias_del_agente[t_omv].costo_om_variable)
                end
            end
        end
        if !isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
            for t_omv_alm in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE
                 for (e_omv,h_omv) in T_global
                    add_to_expression!(costo_om_var_expr_terms, descarga_ag[t_omv_alm,(a_omv,m_omv),(e_omv,h_omv)] * tecnologias_del_agente[t_omv_alm].costo_om_variable)
                end
            end
        end
    end
    coste_omv_ag_total = @expression(model, costo_om_var_expr_terms) # Este es el costo O&M variable total para el NPV
    # Para costo anualizado, se debe hacer por año dentro del bucle de ANIOS si es necesario para otros cálculos anuales.
    # O se puede definir una expresión indexada por año:
    coste_omv_ag_anual = @expression(model, [p in PERIODOS],
        sum(produccion_ag[t,p,(e,h)]*tecnologias_del_agente[t].costo_om_variable for t in TECNOLOGIAS_OPERABLES_AGENTE if !(t in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) for (e,h) in T_global; init=0.0) +
        sum(descarga_ag[t,p,(e,h)]*tecnologias_del_agente[t].costo_om_variable for t in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE for (e,h) in T_global; init=0.0)
    )

    coste_compra_carga_ag = @expression(model, [p in PERIODOS],
    isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) ? 0.0 :
    sum(
        carga_ag[tec, p, (e,h)] *
        get(precios_mercado_esperados,
        (p[1], tecnologias_del_agente[tec].vector, e, h), 0.0)
        for tec in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
            (e,h) in T_global ; init = 0.0)
    )

    coste_co2_ag = @expression(model, [p in PERIODOS],
        (sum(produccion_ag[t,p,(e,h)]*tecnologias_del_agente[t].emisiones_unitarias
            for t in TECNOLOGIAS_OPERABLES_AGENTE if !(t in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE)
            for (e,h) in T_global; init=0.0)
        ) * get_precio_co2_val(p[1])
    )

    coste_total_anual_ag = @expression(model, [p in PERIODOS], coste_inv_ag[p] + coste_omf_ag[p] + coste_omv_ag_anual[p] +
        coste_co2_ag[p] + coste_compra_carga_ag[p])


    ingreso_total_anual_ag = @expression(model, [p in PERIODOS],
        sum(
            (produccion_ag[tec_ing, p, (est_ing,hor_ing)] +
             (tec_ing in TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE ? descarga_ag[tec_ing, p, (est_ing,hor_ing)] : 0.0)
            ) *
            get(precios_mercado_esperados, (p[1], tecnologias_del_agente[tec_ing].vector, est_ing, hor_ing), 0.0)
            for tec_ing in TECNOLOGIAS_OPERABLES_AGENTE, (est_ing,hor_ing) in T_global; init=0.0
        )
    )



            #    REVISAR: ME LAS HE CARGADO (LAS DOS PRIMERAS) PORQUE CREO QUE ME HACEN INFEASIBLE EL MODELO A PARTIR DEL MES 2
    # ─── Restricciones financieras ───────────────────────────────────────
    @constraint(model, [p in PERIODOS], deuda_total_ag[p]  <= MAX_RATIO_DEUDA * patrimonio_neto_ag[p])
    @constraint(model, [p in PERIODOS], cuota_total_ag[p]  <= DSCR_MAX * ingreso_total_anual_ag[p])
    # @constraint(model, [(a,m) in PERIODOS], cap_total_mw_ag[(a,m)] >= CAP_MIN_MW)
    ##################################################################

    aniobase = first(PERIODOS)[1]  

    t0 = time()
    aniobase = first(PERIODOS)[1]

    beneficio_desc = @expression(model, [p in PERIODOS],
        sum((ingreso_total_anual_ag[p] - coste_total_anual_ag[p]) /
            (1 + tasa_descuento)^(p[1] - aniobase)  for p in PERIODOS))

    penalizacion_desc = @expression(model, [p in PERIODOS],
        PENAL_EXCESO * (
            sum(exceso_cap[t,p]  / (1 + tasa_descuento)^(p[1] - aniobase)
                for t in TECNOLOGIAS_OPERABLES_AGENTE, p in PERIODOS) +
            sum(exceso_cash[p]   / (1 + tasa_descuento)^(p[1] - aniobase)
                for p in PERIODOS)
        ))





    # Límite de inversión basado en cash inicial del agente (ejemplo muy simple)
    presupuesto_inversion_total_agente = agente_obj.cash # Por ejemplo, puede invertir hasta 2 veces su cash inicial en todo el horizonte
    @constraint(model, sum(coste_inv_ag[a_inv_lim] for a_inv_lim in PERIODOS) <= presupuesto_inversion_total_agente)
    # REVISAR, PODRÍA SER INTERESANTE INTRODUCIRLO

    beneficio_obj_ag = @expression(model, [p in PERIODOS],
        (beneficio_desc - penalizacion_desc) / 1e6
    )
    
    return (
        id_agente        = id_agente,
        beneficio        = beneficio_obj_ag,
        produccion       = produccion_ag,
        cap_total        = cap_total_ag,
        cap_nueva        = cap_nueva_ag,
        carga            = isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) ? nothing : carga_ag,
        descarga         = isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) ? nothing : descarga_ag,
        nivel            = isempty(TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE) ? nothing : nivel_ag,
        tecs_op          = TECNOLOGIAS_OPERABLES_AGENTE,
        tecs_alm         = TECNOLOGIAS_ALMACENAMIENTO_DEL_AGENTE,
        PERIODOS         = PERIODOS,
        T_global         = T_global,
    )
end

include("PlanificacionMes.jl")


include("PoliticaVerde.jl"); using .PoliticaVerde

#include("MercadoDiario.jl")

    
function resolver_despacho!(model::JuMP.Model)
    JuMP.optimize!(model)
end




# --- helpers para que otros módulos accedan rápido a variables del modelo ---
var_cap_nueva(m) = m[:cap_nueva]
var_cap_total(m) = m[:cap_total]
var_gen(m)       = m[:gen]
emisiones_unitarias(m) = m[:emisiones_unitarias]      # Param. registrado al crear el MILP

function update_mes!(model, demanda_mes::DataFrame,
    exp_eolica::Dict{Int,Float64})
    # Sobrescribe parámetros JuMP `DEM[hora]` y `EOL[h]`
    for (i, row) in enumerate(eachrow(demanda_mes))
        t = i - 1  # timestep que va de 0 a horizon_months*24-1
        JuMP.fix(model[:DEM][t], row.demanda_mw; force = true)
        JuMP.fix(model[:EOL][t], get(exp_eolica, t, 1.0); force = true)
    end
end

set_precio_co2!(model, precio) =
    JuMP.fix(model[:pco2], precio; force = true)
function postprocess_mes_capex!(model, agentes)
    for ag in keys(agentes), tecnologia in keys(agentes[ag].capacidades)
        extra = JuMP.value(var_cap_nueva(model)[ag, tecnologia, 1])  # m=1 es el mes actual
        agentes[ag].capacidades[tecnologia] += extra
        agentes[ag].debt += extra *
            (1 - EQ_SHARE) *
            tecnologias_del_agente[tecnologia].costo_inicial
    end
end

function set_cv_offerta!(model, ag, tecnologia, nuevo_cv)
    JuMP.fix(model[:cv_offerta][ag, tecnologia], nuevo_cv; force = true)
end

function update_dia!(model, demanda_df::DataFrame, perfiles::Dict)
    for h in 0:23
        JuMP.fix(model[:DEM][h], demanda_df.demanda_mw[h]; force = true)
        for (tecnologia, fac) in perfiles
            JuMP.fix(model[:PERF][tecnologia, h], fac[h]; force = true)
        end
    end
end

function postprocess_dia!(model, agentes, precio_co2, subv_unit)
    for ag in keys(agentes)
        prof = JuMP.value(model[:beneficio][ag])
        emis = JuMP.value(model[:emis][ag])
        agentes[ag].profit    += prof
        agentes[ag].emissions += emis
    end
end

function acumular_gen_verde!(acum::Dict{String,Float64}, 
    agentes::Dict{Int,Main.Agentes.EnergyAgent}, 
    fecha::Date)
    # Acumular generación verde del día
    for ag in values(agentes)
        if hasproperty(ag, :gen_diaria) && haskey(ag.gen_diaria, fecha)
            for (tecnologia, gen24h) in ag.gen_diaria[fecha]
                if haskey(agentes[ag.id].var_costs, tecnologia) && 
                agentes[ag.id].var_costs[tecnologia] ≈ 0.0  # criterio "verde"
                acum[tecnologia] = get(acum, tecnologia, 0.0) + gen24h
                end
            end
        end
    end
end



"""
    construir_resultados_agente(
        id_agente,
        tecs_op,
        tecs_alm,
        periodos,
        T_global,
        cap_total,
        cap_nueva,
        produccion,
        carga,
        descarga,
        nivel,
    ) -> (cap_df, prod_df, alm_df)

Devuelve tres `DataFrame`s con la capacidad, la producción y (si procede) el
almacenamiento del agente, construidos de forma vectorizada y sin `push!`.
"""
function construir_resultados_agente(
    id_agente::Int,
    tecs_op::Vector{String},
    tecs_alm::Vector{String},
    periodos,
    T_global,
    cap_total,
    cap_nueva,
    produccion,
    carga,
    descarga,
    nivel,
)
    # ——— Capacidad ————————————————————————————————————————————
    cap_df = DataFrame(
        agent_id  = fill(id_agente, length(tecs_op)*length(periodos)),
        tecnologia = [tec for tec in tecs_op for _ in periodos],
        anio       = [a  for _   in tecs_op for (a,_) in periodos],
        mes        = [m  for _   in tecs_op for (_,m) in periodos],
        cap_total  = [value(cap_total[tec,p])  for tec in tecs_op for p in periodos],
        cap_nueva  = [value(cap_nueva[tec,p])  for tec in tecs_op for p in periodos],
    ; copycols = false)

    # ——— Producción ——————————————————————————————————————————
    prod_df = DataFrame(
        agent_id  = fill(id_agente, length(tecs_op)*length(periodos)*length(T_global)),
        tecnologia = [tec for tec in tecs_op for _ in periodos for _ in T_global],
        anio       = [a   for _ in tecs_op for (a,_) in periodos for _ in T_global],
        mes        = [m   for _ in tecs_op for (_,m) in periodos for _ in T_global],
        estacion   = [e   for _ in tecs_op for _ in periodos for (e,_) in T_global],
        hora       = [h   for _ in tecs_op for _ in periodos for (_,h) in T_global],
        valor      = [value(produccion[tec,p,(e,h)])
                      for tec in tecs_op for p in periodos for (e,h) in T_global],
    ; copycols = false)

    # ——— Almacenamiento (si hay) ——————————————
    alm_df = DataFrame()
    if !isempty(tecs_alm)
        alm_df = DataFrame(
            agent_id  = fill(id_agente, length(tecs_alm)*length(periodos)*length(T_global)),
            tecnologia = [tec for tec in tecs_alm for _ in periodos for _ in T_global],
            anio       = [a   for _ in tecs_alm for (a,_) in periodos for _ in T_global],
            mes        = [m   for _ in tecs_alm for (_,m) in periodos for _ in T_global],
            estacion   = [e   for _ in tecs_alm for _ in periodos for (e,_) in T_global],
            hora       = [h   for _ in tecs_alm for _ in periodos for (_,h) in T_global],
            carga      = [value(carga[tec,p,(e,h)])
                          for tec in tecs_alm for p in periodos for (e,h) in T_global],
            descarga   = [value(descarga[tec,p,(e,h)])
                          for tec in tecs_alm for p in periodos for (e,h) in T_global],
            nivel      = [value(nivel[tec,p,(e,h)])
                          for tec in tecs_alm for p in periodos for (e,h) in T_global],
        ; copycols = false)
    end

    return cap_df, prod_df, alm_df
end



export build_model_capex_base   # ⬅️ añade a la lista de export

"""
build_model_capex_base(agentes, tecnologias; horizon_months=60) -> JuMP.Model

Genera un modelo JuMP mínimo con las variables y parámetros que
`PlanificacionMes` necesita:  
  • cap_nueva, cap_total, gen  
  • parámetros DEM, EOL y pco2

No incluye restricciones, sólo sirve como *plantilla* ligera que luego se
clona y se completa cada mes.
"""
function build_model_capex_base(
        agentes::Dict{Int,Agentes.EnergyAgent},
        tecnologias::Dict{String,Tecnologia};
        horizon_months::Int = 60)

    AG   = collect(keys(agentes))
    all_tecs = union(keys(tecnologias),                    # TODAS las tecnologías definidas
                 (k for ag in values(agentes)   # + las que ya poseen los agentes
                      for k in keys(ag.capacidades)))
    TEC  = collect(all_tecs)
    MES  = 1:horizon_months
    HORAS = 0:(24*horizon_months-1)

    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "threads", Threads.nthreads())
    set_optimizer_attribute(model, "parallel", "on")
    set_optimizer_attribute(model, "solver", "simplex")

    # ── Variables de capacidad ───────────────────────────────────────────────
    @variable(model, cap_nueva[ag in AG, t in TEC, m in MES] ≥ 0)
    @variable(model, cap_total[ag in AG, t in TEC, m in MES] ≥ 0)
    # ── Variable de generación agregada (se detallará después) ───────────────
    @variable(model, gen[ag in AG, t in TEC, h in HORAS] ≥ 0)
    # ── Parámetros que se actualizan cada mes ────────────────────────────────
    @variable(model, DEM[h in HORAS]  >= 0,  start = 0.0)
    @variable(model, EOL[h in HORAS]  >= 0,  start = 1.0)
    @variable(model, pco2 >= 0, start = 0.0)

    return model
end



# ──────────────────────────────────────────────────────────────────────────────
#  FUNCIONES FINANCIERAS CONSERVADORAS (heredadas de PlanificacionMes)
# ──────────────────────────────────────────────────────────────────────────────

"""
    puede_financiar_inversion(agente, tecnologias, tecnologia, mw)

Devuelve `true` si el agente dispone de caja y ratio de endeudamiento  
aceptables tras realizar la inversión propuesta.
"""
function puede_financiar_inversion(agente::Agentes.EnergyAgent,
                                   tecnologias::Dict{String,Tecnologias.Tecnologia},
                                   tecnologia::String,
                                   mw::Float64)


    tec         = tecnologias[tecnologia]
    capex       = mw * tec.costo_inicial
    equity_need = capex * FinanzasParametros.EQ_SHARE

    patrimonio  = Agentes.calcular_patrimonio_neto(agente)
    deuda_fut   = agente.debt + (capex - equity_need)
    ratio_fut   = patrimonio > 0 ? deuda_fut / patrimonio : Inf

    return (agente.cash - equity_need) > CASH_MIN_OPER && ratio_fut < RATIO_DEUDA_MAX
end


end # module OptimizacionSistema