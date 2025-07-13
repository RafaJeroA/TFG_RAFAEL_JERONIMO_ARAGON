# --- START OF CORRECTED FILE Datos.jl ---
module Datos

using DataFrames, CSV, Unicode
using Dates

# Export necessary functions and types
export cargar_y_preparar_datos_base, demanda,
       normalizar_cabeceras!,
       dividir_demanda_por_sector, 
       aplicar_politicas_industria!, aplicar_politicas_almacenamiento!, perfiles_horarios_estaciones_hora, PARAMS



# (No need to export Tecnologia if defined in Tecnologias.jl)

const TEC_NOTIFICADAS = Set{String}()



const PARAMS = Dict{Symbol,Any}()


#try
#    df = CSV.read("Datos/param_simulacion.csv", DataFrame, header=true,
#                  rename=["clave","valor"])
#    for row in eachrow(df)
#        PARAMS[Symbol(row.clave)] = parse(Float64, row.valor)
#    end
#catch e
#    @error "No se pudo leer param_simulacion.csv: $e"
#end

# Cargar precio_base_deficit
try
    df = CSV.read("Datos/precio_base_deficit.csv", DataFrame)
    rename!(df, [:vector,:precio_base_deficit])       # ajusta si ya vienen así
    PARAMS[:precio_base_deficit] =
        Dict(Symbol.(df.vector) .=> Float64.(df.precio_base_deficit))
catch e
        @error "No se pudo leer precio_base_deficit.csv: $e"
end


# ────────────────────────────────────────────────────────────────────
#  ⚡ NUEVO: Cargar parámetros (a_v , b_v) de la curva inversa Cournot
#           desde  Datos/params_cournot.csv
# ────────────────────────────────────────────────────────────────────
try
    df = CSV.read(joinpath("Datos", "params_cournot.csv"), DataFrame)
    rename!(df, [:vector, :a_v, :b_v]; makeunique=true)   # asegura nombres

    # Guardamos en PARAMS los dos diccionarios que exige EquilibrioMultivector
    PARAMS[:dem_inv_a] = Dict(Symbol.(df.vector) .=> Float64.(df.a_v))
    PARAMS[:dem_inv_b] = Dict(Symbol.(df.vector) .=> Float64.(df.b_v))
    @info "Curvas inversas de demanda cargadas: $(length(df.vector)) vectores."
catch e
    @error "No se pudo leer params_cournot.csv: $e"
end



# Columnas requeridas tras normalizar

# (define req_dem_cols arriba del archivo, si no lo tienes:)
#const REQ_DEM_COLS = [:anio, :sector, :subsector, :tecnologia,
#    :demanda_mwh, :dia_del_anio, :vector]
# Verificar cabeceras calendario (asumiendo originales como "anio", "dia_del_anio", "estacion")
#const req_cal = [:anio, :dia_del_anio, :estacion] # Asumiendo que se normalizan a esto
# Verificar perfiles (asumiendo "tecnologia", "estacion", "hora", "factor")
#const req_perf = [:tecnologia, :estacion, :hora, :factor] # Asumiendo que se normalizan a esto

"""
    generar_calendario(anio_inicio, anio_fin;
                       plantilla = "Datos/calendario_estaciones_base.csv") → DataFrame

Devuelve un DataFrame con las fechas desde `anio_inicio` hasta `anio_fin`.
Si un año no es bisiesto, omite el 29-febrero.
Añade una columna `fecha::Date` ya construida.
"""
function generar_calendario(anio_inicio::Int, anio_fin::Int;
    plantilla::String = joinpath("Datos", "calendario_estaciones_base.csv"))
    cal_base = CSV.read(plantilla, DataFrame)

    # Verificamos que solo haya las columnas esperadas
    # 1) Convertimos los nombres de columna a símbolos en minúsculas
    rename!(cal_base,
        [Symbol(lowercase(string(n))) for n in names(cal_base)])

    # 2) Comprobamos que tengamos exactamente :mes, :dia y :estacion
    let req = [:mes, :dia, :estacion], have = Symbol.(names(cal_base))
        faltan = setdiff(req, have)
        if !isempty(faltan)
          @info "Columnas detectadas en calendario: $have"
          error("Faltan columnas en plantilla de calendario: $(faltan)")
        end
    end

    # Armamos el calendario completo
    cal_out = DataFrame(anio = Int[], mes = Int[], dia = Int[], estacion = String[])
    for anio in anio_inicio:anio_fin
        for row in eachrow(cal_base)
            m, d = row.mes, row.dia
            # Saltamos 29-F en años no bisiestos Y verificamos fecha válida
            if (m == 2 && d == 29 && !isleapyear(anio)) || !checkbounds(Bool, 1:Dates.daysinmonth(anio, m), d)
                continue
            end
            push!(cal_out, (anio, m, d, row.estacion))
        end
    end

    # Construimos la columna Date
    cal_out.fecha = Date.(cal_out.anio, cal_out.mes, cal_out.dia)
    # 👉 Añadir día del año para compatibilidad con el resto del flujo
    cal_out.dia_del_anio = dayofyear.(cal_out.fecha)
    return cal_out
end


"""
Normaliza un nombre de columna a un Symbol en minúsculas y sin caracteres especiales.
Ejemplos:
"Anio" -> :anio
"Dia del Anio" -> :dia_del_anio
"Demanda mWh" -> :demanda_mwh
"Precio co2 (€/t)" -> :precio_co2_eur_t
"Anio" -> :anio
"Dia_del_anio" -> :dia_del_anio
"""
function normalizar_nombre_columna(nombre::AbstractString)
    # 1. Normalización Unicode NFD y quitar diacríticos (acentos, etc.)
    limpio = Unicode.normalize(String(nombre), :NFD)
    limpio = filter(c -> !('\u0300' <= c <= '\u036F'), limpio) # Quitar marcas diacríticas

    # 2. Convertir a minúsculas y reemplazar caracteres no deseados por '_'
    # \p{L} = Letras, \p{N} = Números
    limpio = lowercase(replace(limpio, r"[^\p{L}\p{N}_]" => "_"))

    # 3. Colapsar múltiples underscores seguidos a uno solo
    limpio = replace(limpio, r"__+" => "_")

    # 4. Quitar underscores al principio o al final
    limpio = strip(limpio, '_')

    # 5. Manejar caso especial donde todo se reemplaza (nombre original solo símbolos)
    if isempty(limpio)
        return Symbol("_col_", hash(nombre))
    end

    # 6. Asegurarse que no empiece con número (prefijo '_')
    if !isempty(limpio) && isdigit(first(limpio))
        limpio = "_" * limpio
    end

    # 7. Convertir a Symbol
    return Symbol(limpio)
end


"""
Aplica la normalización de nombres a todas las columnas de un DataFrame IN PLACE.
"""
function normalizar_cabeceras!(df::DataFrame)
    nuevos_nombres = [normalizar_nombre_columna(n) for n in names(df)]
    if length(unique(nuevos_nombres)) != length(nuevos_nombres)
        dups = filter(x -> count(==(x), nuevos_nombres) > 1, unique(nuevos_nombres))
        original_names = names(df)
        mapping = Dict(zip(original_names, nuevos_nombres))
        @warn "Nombres de columna duplicados después de normalización: $dups. Originales -> Normalizados: $mapping. Se usará 'makeunique=true'."
    end
    try
        rename!(df, nuevos_nombres, makeunique=true)
    catch e
        println("Error crítico normalizando cabeceras: $e")
        println("Nombres originales: $(names(df))")
        println("Nombres normalizados intentados: $nuevos_nombres")
        rethrow(e)
    end
end





"""
Carga datos base y genera demanda HORARIA a partir de demanda anual, perfiles estacionales y horarios.
SEPARA perfiles de demanda (para distribución horaria) y perfiles de generación (para disponibilidad).
"""
function cargar_y_preparar_datos_base(; anio_inicio = 2024, anio_fin = 2050)
    println("⏳ Cargando datos base de demanda (anual, estacionales y horarios)...")
    
    # ═══════════════════════════════════════════════════════════════════
    # 1. DATOS BASE Y CALENDARIO
    # ═══════════════════════════════════════════════════════════════════
    df_base = CSV.read(joinpath("Datos","demanda_base.csv"), DataFrame)
    normalizar_cabeceras!(df_base)
    calendario_df = generar_calendario(anio_inicio, anio_fin)

    df_cal = calendario_df
    df_cal.dia_del_anio = dayofyear.(df_cal.fecha)
    df_days = combine(groupby(df_cal, [:anio, :estacion]),
                      :dia_del_anio => length => :n_dias)

    # ═══════════════════════════════════════════════════════════════════
    # 2. PERFILES ESTACIONALES (solo para demanda)
    # ═══════════════════════════════════════════════════════════════════
    df_est = CSV.read(joinpath("Datos","perfiles_estacionales_demanda.csv"), DataFrame)
    normalizar_cabeceras!(df_est)

    # ═══════════════════════════════════════════════════════════════════
    # 3. PERFILES HORARIOS DE DEMANDA (para distribución de demanda anual)
    # ═══════════════════════════════════════════════════════════════════
    hor_demanda_file = joinpath("Datos","perfiles_demanda.csv")
    df_hor_demanda = DataFrame(estacion=String[], hora=Int[])
    try
        if !isfile(hor_demanda_file) || filesize(hor_demanda_file) == 0
            throw(ErrorException("'perfiles_demanda.csv' inexistente o vacío"))
        end
        df_hor_demanda = CSV.read(hor_demanda_file, DataFrame)
        normalizar_cabeceras!(df_hor_demanda)
        if :hour in names(df_hor_demanda)
            rename!(df_hor_demanda, :hour => :hora)
        end
        @info "Perfiles horarios de DEMANDA cargados desde $(basename(hor_demanda_file))"
    catch e
        @warn "No se pudo cargar perfiles_demanda.csv: $e. Generando perfil plano para demanda." 
    end
    
    # Completar vectores de demanda faltantes con perfil plano
    expected_demand_vectors = intersect(setdiff(names(df_base), [:anio]),
                                       setdiff(names(df_est), [:estacion]))
    faltantes_demanda = setdiff(expected_demand_vectors, names(df_hor_demanda))
    for v in faltantes_demanda
        df_hor_demanda[!, v] = fill(1/24, nrow(df_hor_demanda))
    end

    # ═══════════════════════════════════════════════════════════════════
    # 4. PERFILES HORARIOS DE GENERACIÓN (para disponibilidad de tecnologías)
    # ═══════════════════════════════════════════════════════════════════
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

    # ═══════════════════════════════════════════════════════════════════
    # 5. PROCESAR DEMANDA HORARIA (usando perfiles de demanda)
    # ═══════════════════════════════════════════════════════════════════
    df_base_long = stack(df_base, Not(:anio), variable_name=:vector, value_name=:demanda_anual)
    df_est_long  = stack(df_est,  Not(:estacion), variable_name=:vector, value_name=:factor_est)
    df_hor_demanda_long = stack(df_hor_demanda, Not([:estacion,:hora]), variable_name=:vector, value_name=:factor_hor)

    # Demanda diaria por estación
    df_base_est = innerjoin(df_base_long, df_est_long, on=:vector)
    df_base_est = innerjoin(df_base_est, df_days, on=[:anio,:estacion])
    transform!(df_base_est, [:demanda_anual,:factor_est,:n_dias] => ByRow((d,f,n)-> d*f/n) => :demanda_diaria)
    

    # Distribuir demanda horaria según perfiles horarios DE DEMANDA
    if nrow(df_hor_demanda_long) == 0
        println("⚠️ No hay perfiles horarios de demanda definidos. Generando perfil plano uniforme.")
        filas = []
        sizehint!(filas, nrow(df_base_est)*24)
        for row in eachrow(df_base_est), h in 0:23
            push!(filas, (
                anio     = Int(row.anio),
                estacion = String(row.estacion),
                hora     = h,
                vector   = String(row.vector),
                valor    = round(row.demanda_diaria/24, digits=6)
            ))
        end
        df_horaria = DataFrame(filas)
    else
        df_horaria = innerjoin(df_base_est, df_hor_demanda_long, on=[:estacion,:vector])
        transform!(df_horaria, [:demanda_diaria,:factor_hor] => ByRow(*) => :valor)
        df_mes = unique(select(calendario_df, [:anio, :estacion, :mes]))
        df_horaria = leftjoin(df_horaria, df_mes, on = [:anio, :estacion])
        
        # Reordenar para coherencia
        select!(df_horaria,
                [:anio, :mes, :estacion, :hora, :vector, :valor])
    end
    println("✅ Demanda horaria generada con $(nrow(df_horaria)) filas.")

    # ═══════════════════════════════════════════════════════════════════
    # 6. CREAR DICCIONARIO DE PERFILES DE GENERACIÓN 
    # ═══════════════════════════════════════════════════════════════════
    perfiles_gen = Dict{Tuple{String,String},Dict{Int,Float64}}()
    
    # ✅ USAR df_hor_generacion (tecnologías productoras), NO df_hor_demanda_long
    for row in eachrow(df_hor_generacion)
        key = (String(row.tecnologia), String(row.estacion))  # ← TECNOLOGÍA, no vector
        if !haskey(perfiles_gen, key)
            perfiles_gen[key] = Dict{Int,Float64}()
        end
        perfiles_gen[key][Int(row.hora)] = Float64(row.factor)
    end
    
    println("✅ Perfiles de generación procesados: $(length(perfiles_gen)) combinaciones tecnología-estación.")

    # ═══════════════════════════════════════════════════════════════════
    # 7. OTROS DATOS
    # ═══════════════════════════════════════════════════════════════════
    pot_df   = CSV.read("Datos/potencial_renovable.csv", DataFrame); normalizar_cabeceras!(pot_df)
    delay_df = CSV.read("Datos/delay_construccion.csv", DataFrame); normalizar_cabeceras!(delay_df)

    pot_dict   = Dict(row.tecnologia => row.potencial_max_mw for row in eachrow(pot_df))
    delay_dict = Dict(row.tecnologia => row.delay             for row in eachrow(delay_df))
    
    cal_estaciones = Dict(r.fecha => r.estacion for r in eachrow(calendario_df))
    exp_eolica = Dict{Date,Float64}()

    return (
        demanda_df          = df_horaria,              # ← Demanda distribuida horariamente
        potencial_renovable = pot_dict,
        delay_construccion  = delay_dict,
        cal_estaciones      = cal_estaciones,
        perfiles_gen        = perfiles_gen,            # ← Dict para optimización (tecnología,estación) → horas
        perfiles_generacion = df_hor_generacion,      
        exp_eolica          = exp_eolica,
    )
end


# 📥 Cargar demanda HORARIA
const demanda = try
                    cargar_y_preparar_datos_base()
                catch e
                    @error "FALLO CRÍTICO AL CARGAR datosNT.demanda_df: $e"
                    # Re-throw original exception para mostrar stacktrace
                    rethrow(e)
                end

# --- Funciones para Segmentar y Aplicar Políticas ---
# (Usan nombres normalizados)

function dividir_demanda_por_sector(df::DataFrame)
    if isempty(df)
        @warn "Intentando dividir un DataFrame de demanda vacío."
        return (transporte=DataFrame(), industria=DataFrame(), almacenamiento=DataFrame(), otros=DataFrame())
    end
    if :sector ∉ Symbol.(names(df))
        @warn "La columna ':sector' no existe. Asignando todo a 'otros'."
        return (transporte=DataFrame(), industria=DataFrame(), almacenamiento=DataFrame(), otros=df)
    end

    sectores_str = String.(coalesce.(df.sector, "")) # Manejar missings
    filter_transporte = lowercase.(sectores_str) .== "transporte"
    filter_industria = lowercase.(sectores_str) .== "industria"
    filter_almacenamiento = lowercase.(sectores_str) .== "almacenamiento"

    df_transporte = df[filter_transporte, :]
    df_industria = df[filter_industria, :]
    df_almacenamiento = df[filter_almacenamiento, :]
    df_otros = df[.!filter_transporte .& .!filter_industria .& .!filter_almacenamiento, :]

    return (transporte = df_transporte, industria = df_industria, almacenamiento = df_almacenamiento, otros = df_otros)
end




function fechas_del_mes(a::Int, m::Int)::Vector{Date}
    inicio = Date(a, m, 1)
    return collect(inicio:Dates.lastdayofmonth(inicio))
end

function demanda_dia(datos, f::Date)
    # Obtener la estación de la fecha usando el calendario
    estacion_fecha = get(datos.cal_estaciones, f, nothing)
    if isnothing(estacion_fecha)
        @warn "No se encontró estación para la fecha $f en el calendario"
        return DataFrame()  # Retornar DataFrame vacío
    end
    
    # Filtrar demanda por año y estación
    filter(r -> r.anio == year(f) && r.estacion == estacion_fecha,
        datos.demanda)
end

function perfiles_dia(datos, f::Date)
    est = get(datos.cal_estaciones, f, nothing)
    if est === nothing
        @warn "No hay estación para $f"
        return Dict()
    end
    return Dict(
      vec => horario
      for ((vec, estacion), horario) in datos.perfiles_gen
      if estacion == est
    )
end



function perfiles_horarios_estaciones_hora()
    # Devuelve lista de tuplas (estación, hora) únicas
    return unique([(String(r.estacion), Int(r.hora)) 
        for r in eachrow(demanda.demanda_df)]) #PUEDE QUE DEMANDA NO SEA DATOS.DEMANDA, CUIDADO
end

const Bundle = NamedTuple


"Versión rápida para todo el día: devuelve Dict{Int,Dict{String,Float64}} hora→(vector→mw)"
function demanda_dia_todos_vectores(datos::Bundle, fecha::Date)
    # obtenemos la estación (verano/invierno) para esa fecha
    est = get(datos.cal_estaciones, fecha, nothing)
    if est === nothing
        @warn "No hay estación definida para la fecha $fecha"
        return Dict{Int,Dict{String,Float64}}()
    end
    # filtramos el perfil horario por esa estación
    df = @view datos.demanda[datos.demanda.estacion .== est, :]


    out = Dict{Int,Dict{String,Float64}}()
    for g in groupby(df, [:hora, :vector])
        h = Int(first(g.hora))
        v = first(g.vector)
        out[h] = get(out, h, Dict{String,Float64}())
        out[h][v] = first(g.valor)
    end
    return out
end                           # >>> NUEVO

end  # module Datos
