# =========================================
# main.jl  –  Entrada única del modelo
# =========================================
# Ejecuta escenarios de simulación energética integrando
# optimización MILP, ABM, análisis de redes, dinámica y juegos.
#
# USO:
#   julia main.jl             # Ejecución secuencial
#   julia main.jl --parallel  # Ejecución paralela usando workers
# -----------------------------------------

# ────────────── Opciones CLI ─────────────
# Ruta al DLL que acabas de compilar
ENV["HIGHS_LIBRARY_PATH"] = raw"C:\Users\rafae\Desktop\TFG\Modelo_2.5.1\HiGHS\build\Release\bin\highs.dll"
ENV["CUDA_HOME"]         = raw"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.9"
ENV["HIGHS_CUDA_DEVICE"]  = "0"          # GPU que quieres usar
ENV["OMP_NUM_THREADS"]    = "27"         # hilos para la parte CPU
const PARALLEL = true  # Always parallel to ensure worker functions are defined

using Logging
#global_logger(SimpleLogger(stderr, Logging.Info))
# ────────────── Dependencias Principales ───────
using Distributed, Dates, CSV, DataFrames, Plots, StatsPlots # Plots/StatsPlots para Resultados
using JuMP, HiGHS # Para la optimización (HiGHS es el solver)
using DifferentialEquations # Para DinamicaHidrogeno
try
    using Agents # Para ABM (puede fallar en Windows con ReadOnlyMemoryError)
catch e
    @warn "Fallo cargando paquete Agents: $e. Se omitirá ABM externo."
end
using MathOptInterface # Para MOI status
const MOI = MathOptInterface
# using GameTheory # Descomentar si se usa activamente JuegoEstrategico y está instalado
using XLSX
global datosNT

using Dates


# Declaración explícita de variables globales
global resultados_completos = Vector{Any}()
global resultados_ok        = 0
global errores_post         = 0

const n_workers = 1
const HIGHS_THREADS = Threads.nthreads()/n_workers

#global  resultados_completos = Vector{Any}(), resultados_ok=0, errores_post=0


# LIMITACIÓN PARA QUE NO EXPLOTE EL ORDENADOR
"""
Devuelve el nº «seguro» de workers ≈ 75 % de los núcleos lógicos,
dejando al menos 1 hilo libre para el SO y nunca menos de 1 worker.

MI ORDENADOR NO PUEDE CORRERLO, POR LO QUE NO LO USO, PERO ES POTENCIALMETE ÚTIL
"""
function workers_75_pct()
    ncpu = Sys.CPU_THREADS              # núcleos lógicos
    n    = max(1, floor(Int, 0.75 * ncpu))
    return n
end

println("Julia Threads: ", Threads.nthreads(), " | CPU cores: ", Sys.CPU_THREADS)


function normalizar_cabeceras_worker!(df::DataFrame)
    nuevos_nombres = [normalizar_nombre_columna_worker(n) for n in names(df)]
     try rename!(df, nuevos_nombres, makeunique=true) catch e
         println("Error worker normalizando cabeceras: $e"); println("Nombres actuales: $(names(df))"); println("Intentados: $nuevos_nombres");
     end
end



# ────────────── Módulos del Proyecto (ORDENADOS POR DEPENDENCIA) ─────────
println("🔄 Cargando módulos del proyecto...")
try
    # 1. Módulos base sin dependencias internas complejas o que definen tipos usados por otros
    include("Escenarios.jl"); using .Escenarios          # Define Escenario struct
    include("Utils.jl"); using .Utils                   # Usa Escenario
    include("Tecnologias.jl"); using .Tecnologias        # Define Tecnologia struct
    include("InversionesPendientes.jl"); using .InversionesPendientes

    # 2. Módulos con dependencias de los anteriores
    include("Datos.jl"); using .Datos
    global datosNT = Datos.cargar_y_preparar_datos_base()
    if isempty(datosNT.demanda_df)
        @error "El DataFrame datosNT.demanda_df está vacío tras cargar Datos.jl. Verifica CSV y cargar_y_preparar_datos_base()."
        exit(1)
    end
    println("   datosNT.demanda_df cargado con $(nrow(datosNT.demanda_df)) filas")

    import .Datos: normalizar_cabeceras!




    
    # 3. Módulo principal de optimización (depende de varios anteriores)
    include("FinanzasParametros.jl"); using .FinanzasParametros
    include("Agentes.jl"); using .Agentes
    include("PoliticaVerde.jl"); using .PoliticaVerde
    include("OptimizacionSistema.jl"); using .OptimizacionSistema # Usa Tecnologias, Escenarios, Datos
    #include("JuegoEstrategico.jl"); using .JuegoEstrategico
    include("PlanificacionMes.jl");   using .PlanificacionMes
    #include("MercadoDiario.jl");   using .MercadoDiario
    include("Resultados.jl");   using .Resultados
    include("SimulacionSecuencial.jl"); using .SimulacionSecuencial


    # include("VisualizacionPlotly.jl") # Descomentar si se prefiere Plotly

catch e
    println("❌ Error fatal al cargar módulos iniciales: $e")
    showerror(stdout, e, catch_backtrace()); println()
    exit(1) # Terminar si los módulos básicos no cargan
end
println("✅ Módulos cargados.")


# ──────────────  Paralelización  ─
if PARALLEL
    n_cores = Sys.CPU_THREADS                    # núcleos lógicos
    n_procs_to_add = n_workers          # Son 3 porque hay 3 escenarios actualmente, CAMBIAR EN UN FUTURO
    @info "Arrancando $n_procs_to_add workers (uno por escenario)…"
    #n_procs_to_add = max(1, floor(Int, 0.75*n_cores))           # ESTO ES DEMASIADO CONSUMO, VOY A HACER QUE SOLAMENTE SE EJECUTEN EL NUMERO DE WORKERS = ESCENARIOS
    #@info "Arrancando $n_procs_to_add workers (≈ 75 % de $n_cores hilos lógicos)…"
    if nprocs() < n_procs_to_add + 1
        try
            addprocs(n_procs_to_add; exeflags=["--project", "--threads", string(HIGHS_THREADS)]) # Asegurar que usen el mismo proyecto
            @info "Modo paralelo: Añadidos $n_procs_to_add workers. Total procesos: $(nprocs())"
        catch e
            @error "Error añadiendo workers: $e. Reintentando con menos workers..."
            try addprocs(max(1, floor(Int, 0.75*n_procs_to_add)); exeflags=["--project", "--threads", string(HIGHS_THREADS)])  ; @info "Workers añadidos: $(nprocs()-1)" 
            catch e
                @error "Error añadiendo workers: $e. Reintentando con menos workers..."
                try addprocs(max(1, n_procs_to_add ÷ 2); exeflags=["--project", "--threads", string(HIGHS_THREADS)])   ; @info "Workers añadidos: $(nprocs()-1)" catch; @warn "No se pudieron añadir workers." end
            end
        end
    else
        @info "Modo paralelo: Ya hay $(nprocs()) procesos activos."
    end

    # Código esencial que debe estar disponible en TODOS los procesos
    println("⏳ Configurando workers paralelos...")
    @everywhere if myid() != 1 begin
        # Sólo los workers (ID ≠ 1) incluyen de nuevo los ficheros:
        include("Escenarios.jl")
        include("Utils.jl")
        include("Tecnologias.jl")
        include("InversionesPendientes.jl")
        include("Datos.jl")
        include("FinanzasParametros.jl")
        include("Agentes.jl")
        include("PoliticaVerde.jl")
        include("OptimizacionSistema.jl")
        include("PlanificacionMes.jl")
        #include("MercadoDiario.jl")
        include("Resultados.jl")
        include("SimulacionSecuencial.jl")
    end

    # Y luego, en todos los procesos (master + workers), traemos los módulos al namespace:
    @everywhere using Main.Escenarios
    @everywhere using Main.Utils
    @everywhere using Main.Tecnologias
    @everywhere using Main.Datos
    @everywhere using Main.Agentes
    @everywhere using Main.PoliticaVerde
    @everywhere using Main.FinanzasParametros
    @everywhere using Main.OptimizacionSistema
    @everywhere using Main.PlanificacionMes
    # @everywhere using Main.MercadoDiario
    @everywhere using Main.Resultados
    @everywhere using Main.SimulacionSecuencial

    # Dependencias básicas adicionales para todos los procesos
    @everywhere using Dates, CSV, DataFrames, Unicode, JuMP
    @everywhere using MathOptInterface # Para MOI status
    @everywhere const MOI = MathOptInterface

    # --- Funciones Helper para Configuración en Workers ---

    # Re-definir normalización aquí para que esté disponible en el worker scope
    function normalizar_nombre_columna_worker(nombre::AbstractString)
         limpio = Unicode.normalize(String(nombre), :NFD)
         limpio = filter(c -> !('\u0300' <= c <= '\u036F'), limpio)
         limpio = lowercase(replace(limpio, r"[^\p{L}\p{N}_]" => "_"))
         limpio = replace(limpio, r"__+" => "_")
         limpio = strip(limpio, '_')
         if isempty(limpio) return Symbol("_col_", hash(nombre)) end
         if !isempty(limpio) && isdigit(first(limpio))
            limpio = "_" * limpio
        end
         return Symbol(limpio)
    end



    # Función para construir el bundle de datos DENTRO del worker
    function construir_data_bundle_worker(agents::Dict{Int,Main.Agentes.EnergyAgent})
        worker_id = myid()
        println("Worker $worker_id: Cargando y preparando datos...")

        # 1. Demanda base horaria (igual que antes)
        demanda_df_horaria = Main.Datos.cargar_y_preparar_datos_base()
            if isempty(demanda_df_horaria) error("Worker $(worker_id): Datos.cargar_y_preparar_datos_base() devolvió un DataFrame vacío.") end

            # Cargar otros archivos necesarios
            perfiles_gen_df = CSV.read(joinpath("Datos", "perfiles_horarios.csv"), DataFrame)
            politicas_df = CSV.read(joinpath("Datos", "politicas.csv"), DataFrame)

            # Normalizar cabeceras de archivos cargados aquí
            normalizar_cabeceras_worker!(perfiles_gen_df)
            normalizar_cabeceras_worker!(politicas_df)

        # --- NUEVO: cap_existente_dict a partir de ``agents`` --------------------
        anio0 = minimum(demanda_df_horaria.anio)               # primer año del horizonte
        cap_existente_dict = Dict{String,Dict{Tuple{Int,Int},Float64}}()
        for ag in values(agents), (tecnologia,cap) in ag.capacidades
            d = get!(cap_existente_dict, tecnologia, Dict{Int,Float64}())
            d[anio0] = get(d, anio0, 0.0) + cap                # acumular mw
        end
        # ------------------------------------------------------------------------

        println("Worker $worker_id: Data bundle construido.")
        return (
            demanda_base_horaria = demanda_df_horaria,
            cap_existente        = cap_existente_dict,
            perfiles_generacion  = perfiles_gen_df,
            politicas            = politicas_df,
            tasa_descuento       = 0.07,
        )
    end


    function construir_tecnologias_bundle_worker()
        worker_id = myid()
        println("Worker $worker_id: Cargando tecnologías...")
        try
             dftec = CSV.read(joinpath("Datos", "tecnologias.csv"), DataFrame)
             normalizar_cabeceras_worker!(dftec)
             tecs = Main.Tecnologias.inicializar_tecnologias!(dftec)
             
             # VERIFICACIÓN CRÍTICA:
            if !haskey(tecs, "coche_gasolina")
                 @error "Worker $worker_id: 'coche_gasolina' NOT found in tecnologias dict!"
                 @error "Worker $worker_id: Available technologies: $(sort(collect(keys(tecs))))"
                 error("Critical technology missing in worker $worker_id")
            else
                 println("Worker $worker_id: ✓ 'coche_gasolina' verified in tecnologias dict")
            end
             
             println("Worker $worker_id: Tecnologías bundle construido with $(length(tecs)) technologies.")
             return tecs
        catch e
             println("❌ Worker $(worker_id): Error construyendo tecnologías bundle: $e")
             showerror(stdout, e, catch_backtrace()); println()
             rethrow(e)
        end
    end

    
end # fin @everywhere
println("✅ Workers configurados.")
end # fin if PARALLEL



# ─────────── Construir bundle de datos globales para competencia ──────────
println("🔄 Preparando datos globales para simulación competitiva…")
# Ya tienes datosNT.demanda_df tras cargar Datos.jl
perfiles_df = CSV.read(joinpath("Datos","perfiles_horarios.csv"), DataFrame)
Datos.normalizar_cabeceras!(perfiles_df)
politicas_df = CSV.read(joinpath("Datos","politicas.csv"), DataFrame)
Datos.normalizar_cabeceras!(politicas_df)

escenarios_definidos = Escenarios.definir_escenarios()

# ───────────────── RELLENAR POLÍTICAS EN LOS ESCENARIOS ─────────────────
function precargar_politicas_en_escenario!(
        esc::Escenarios.Escenario,
        politicas::DataFrame)

    for row in eachrow(politicas)
        anio = Int(row.anio)
        p_co2   = Float64(row.precio_co2)
        subvmax = Float64(row.subv_max_anual)

        # Saltar carga para BAU → todo queda en cero
        if lowercase(esc.nombre) == "bau"
            p_co2   = 0.0
            subvmax = 0.0
        end    

        if lowercase(esc.nombre) == "esfuerzo_x2"
            p_co2   *= 2.0
            subvmax *= 2.0
        elseif lowercase(esc.nombre) == "esfuerzo_mitad"
            p_co2   *= 0.5
            subvmax *= 0.5
        end
        # Precio de CO₂ → los 12 meses del año
        for mes in 1:12
            esc.precio_co2[(anio, mes)] = p_co2
        end
        # Límite anual de subvención
        esc.subv_max_anual[anio] = subvmax
    end
end

for esc in escenarios_definidos
    precargar_politicas_en_escenario!(esc, politicas_df)
end
println("📋   Políticas cargadas en los escenarios: CO₂ y subvención OK.")
# ─────────────────────────────────────────────────────────────────────────



const datos_bundle_main = (
    demanda         = datosNT.demanda_df,
    perfiles_generacion = perfiles_df,
    politicas          = politicas_df,
    tasa_descuento     = 0.07,
    cal_estaciones     = datosNT.cal_estaciones,
    potencial_renovable = datosNT.potencial_renovable,
    delay_construccion  = datosNT.delay_construccion,
    perfiles_gen       = datosNT.perfiles_gen,
    exp_eolica         = datosNT.exp_eolica,
)
println("✅ Data bundle competitivo listo.")




# ─────────── Configuración de Escenarios ─────────
println("📝 Definiendo escenarios...")
#escenarios_definidos = Escenarios.definir_escenarios()
# Crear parámetros para cada escenario
escenarios_para_ejecutar = [
    Dict(
        # Info del escenario
        :nombre_escenario => escena.nombre,
        :escenario_obj    => escena, # Pasar el objeto escenario completo
        :agents_file      => joinpath("Datos","agentes_generacion.csv"),
        :finanzas_iniciales_file => joinpath("Datos","agentes_finanzas.csv"),
        # Parámetros específicos para otros módulos (ajustar según necesidad)
        #:nprod => 10, :ncons => 80, :pasos => 120, # ABM
        #:pay1 => [4.0 1.0; 0.0 2.0], :pay2 => [2.0 0.0; 1.0 3.0]  # JuegoEstrategico
    )
    for escena in escenarios_definidos # Iterar sobre los escenarios definidos
]



    ###################################################
    ###################################################
    #PROVISIONAL PARA HACER LAS PRUEBAS
    escenarios_para_ejecutar = filter(
        par -> par[:nombre_escenario] == "Base",
        escenarios_para_ejecutar
    )
    ###################################################
    ###################################################        





println("-> $(length(escenarios_para_ejecutar)) escenarios listos para ejecutar: $([e[:nombre_escenario] for e in escenarios_para_ejecutar])")


# ─────────── Ejecución de SIMULACIÓN COMPETITIVA (Cournot) ─────────
println("\n🚀 Iniciando simulación competitiva para $(length(escenarios_para_ejecutar)) escenarios...")
tiempo_inicio = now()



# Cargar tecnologías una vez en el proceso principal para pasar a graficar/analizar
tecs = try
    dftec_main = CSV.read(joinpath("Datos", "tecnologias.csv"), DataFrame)
    Datos.normalizar_cabeceras!(dftec_main)
    Tecnologias.inicializar_tecnologias!(dftec_main)
catch e
    @error "Error cargando tecnologías en el proceso principal: $e. Gráficos pueden fallar."
    Dict{String, Main.Tecnologias.Tecnologia}() # Devolver vacío con tipo correcto
 end
# Añadir debug:
#println("🔍 DEBUG: Tecnologías cargadas:")
#for (nombre, tec) in tecs
#    println("  - $nombre")
#end
#println("🔍 Total tecnologías: $(length(tecs))")
#println("🔍 ¿Existe 'coche_gasolina'? $(haskey(tecs, "coche_gasolina"))")



# Función helper para el nuevo modelo (después de línea ~260)
function ejecutar_mercado_diario_mensual(par::Dict)
    println("🏁 Iniciando simulación mercado diario para $(par[:nombre_escenario])...")
    
    # 1) Cargar agentes y finanzas (igual que antes)
    agentes = Agentes.cargar_agentes_generacion(par[:agents_file])
    Agentes.cargar_finanzas_iniciales!(agentes; ruta=par[:finanzas_iniciales_file])
    

    println("🔧 VERIFICACIÓN EN MAIN: Sincronizando var_costs...")
    for (ag_id, agente) in agentes
        tecnologias_faltantes = String[]
        
        for tecnologia in keys(agente.capacidades)
            if !haskey(agente.var_costs, tecnologia)
                push!(tecnologias_faltantes, tecnologia)
                
                # Buscar en tecs con nombre original y normalizado
                if haskey(tecs, tecnologia)
                    agente.var_costs[tecnologia] = tecs[tecnologia].costo_om_variable
                    agente.emissions[tecnologia] = tecs[tecnologia].emisiones_unitarias
                else
                    # Probar con nombre normalizado
                    tecnologia_norm = Utils.normalizar_nombre_tecnologia(tecnologia)
                    if haskey(tecs, tecnologia_norm)
                        agente.var_costs[tecnologia] = tecs[tecnologia_norm].costo_om_variable
                        agente.emissions[tecnologia] = tecs[tecnologia_norm].emisiones_unitarias
                    else
                        agente.var_costs[tecnologia] = 0.0
                        agente.emissions[tecnologia] = 0.0
                        @warn "MAIN: Tecnología '$tecnologia' no encontrada en tecs"
                    end
                end
            end
        end
        
        if !isempty(tecnologias_faltantes)
            println("  ⚠️  MAIN: Agente $(agente.name) sincronizó $(length(tecnologias_faltantes)) tecnologías: $tecnologias_faltantes")
        end
    end

    # 2) Llamar al nuevo módulo SimulacionSecuencial (modificado)
    resultado = SimulacionSecuencial.ejecutar_mercado_real(
        par[:escenario_obj],
        agentes,
        datos_bundle_main,
        tecs
    )
    
    return resultado
end

# Ejecuta secuencialmente todos los escenarios
resultados_completos = [ ejecutar_mercado_diario_mensual(par) for par in escenarios_para_ejecutar ]

tiempo_fin = now()
duracion = canonicalize(round(tiempo_fin - tiempo_inicio, Dates.Second))
println("\n🏁 Simulación competitiva COMPLETADA.")
println("⏱️  Duración total: $duracion")

# ─────────── Procesamiento de Resultados ─────────
println("\n💾 Procesando y guardando resultados...")


path_salidas = "Salidas"
isdir(path_salidas) || mkdir(path_salidas) # Crear directorio si no existe



# Iterar sobre los resultados obtenidos
for (i, res) in enumerate(resultados_completos)
    #if isa(res, String)
    #    nombre = escenarios_para_ejecutar[i][:nombre_escenario]
    #    println("🔄 Procesando traza del escenario '$nombre': $res")
    #    df_traza = CSV.read(res, DataFrame)
    #    Resultados.graficar_traza(res)
    #    continue
    #end

    # Asegurarse de que 'res' es el tipo esperado (NamedTuple)
    if !isa(res, NamedTuple)
        @error "Resultado inesperado para índice $i: $(typeof(res)). Se omite."
        global errores_post += 1
        continue
    end

    esc_nombre = get(res, :escenario, "Desconocido_$i")

    # Verificar si hubo error capturado en el worker
    error_info = get(res, :error_info, nothing)
    if !isnothing(error_info)
        # Ya contamos el error en la fase de ejecución, solo loguear aquí
        @error "Procesando resultado con error capturado en worker para '$esc_nombre'."
        # Mostrar el error capturado si se quiere más detalle aquí también
        # showerror(stdout, error_info.exception, error_info.backtrace); println()
        global errores_post += 1 # Contar errores aquí también para el resumen final
        continue # Saltar al siguiente resultado
    end

    # Verificar estado de la optimización
    opt_status = get(res, :opt_status, MOI.OPTIMIZE_NOT_CALLED)
    resultados_opt = get(res, :opt_resultados, Dict())
    #df_cap = get(resultados_opt, :capacidad, DataFrame())
    #df_prod = get(resultados_opt, :produccion, DataFrame())

    # Ahora tus DataFrames están en res.opt_resultados[:capacidad] y [:produccion]
    df_cap = get(res.opt_resultados, :capacidad, DataFrame())
    df_prod = get(res.opt_resultados, :produccion, DataFrame())



    # Estados considerados exitosos
    estados_exitosos = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
    #push!(resultados_completos, res) # CUIDADO CON ESTA LINEA
   
    # --- Nueva condición: si vengo de Cournot (historial no vacío), lo tratamos como éxito también ---
    es_cournot_ok = !haskey(res, :opt_status) || res.opt_status === nothing && !isempty(get(res, :historial, Dict()))
   
    # Estados con solución (aunque no óptima)
    estados_con_solucion = (MOI.TIME_LIMIT, MOI.ITERATION_LIMIT, MOI.NODE_LIMIT, MOI.SOLUTION_LIMIT, MOI.MEMORY_LIMIT, MOI.OBJECTIVE_LIMIT, MOI.INTERRUPTED)

    if opt_status in estados_exitosos || es_cournot_ok
        println("--- Resultados Óptimos/Aceptables para Escenario: $esc_nombre (Estado: $opt_status) ---")
        opt_obj_val = get(res, :opt_objective, NaN)
        println("   Valor Objetivo: $(round(opt_obj_val, digits=0))")
        global resultados_ok += 1

        if i <= length(escenarios_para_ejecutar)
            esc_obj = escenarios_para_ejecutar[i][:escenario_obj]
        
            # Guardar resultados CSV y Excel
            try
                if !isempty(df_cap) || !isempty(df_prod)
                    # 1) Guardar CSV
                    Resultados.guardar_resultados(df_cap, df_prod, esc_obj)
        
                    # 2) Construir carpeta y nombre de fichero
                    mkpath(path_salidas)
                    # opcional: sanitizar nombre
                    nombre_safe = replace(esc_obj.nombre, r"[\\/:*?\"<>|]" => "_")
                    ruta_xlsx  = joinpath(path_salidas, "resumen_$(nombre_safe).xlsx")
        
                    # 3) Escribir Excel con dos hojas
                    XLSX.openxlsx(ruta_xlsx, mode="w") do xf
                        if !isempty(df_cap)
                            sheet_cap = XLSX.addsheet!(xf, "Capacidad")
                            XLSX.writetable!(sheet_cap, df_cap)
                        end
                        if !isempty(df_prod)
                            sheet_prod = XLSX.addsheet!(xf, "Producción")
                            XLSX.writetable!(sheet_prod, df_prod)
                        end

                        #
                        df_lcoe = get(resultados_opt, :lcoe, DataFrame())
                        if !isempty(df_lcoe)
                            sheet = XLSX.addsheet!(xf, "LCOE")
                            XLSX.writetable!(sheet, df_lcoe)
                        end
                        df_part = get(resultados_opt, :participacion, DataFrame())
                        if !isempty(df_part)
                            sheet = XLSX.addsheet!(xf, "Participacion")
                            XLSX.writetable!(sheet, df_part)
                        end
                        #






                    end
                    println("💾 Excel de resultados guardado en $ruta_xlsx")
                else
                    @warn "No hay datos de capacidad o producción para guardar en escenario $esc_nombre."
                end
            catch e
                @error "Error guardando resultados para $esc_nombre: $e"
                showerror(stdout, e, catch_backtrace()); println()
                global errores_post += 1
            end

            # Graficar resultados
            try
                 Resultados.graficar_resultados(res, esc_obj, tecs)
            catch e
                 @error "Error graficando resultados para $esc_nombre: $e"; showerror(stdout, e, catch_backtrace()); println()
                 global errores_post += 1
            end

            # --- Nueva llamada a guardar_produccion_diaria_ultimo_dia ---
            try
                if !isempty(df_prod)
                    println("   Generando reporte de producción diaria del último día para $esc_nombre...")
                    tecnologias_csv_df = CSV.read(joinpath("Datos", "tecnologias.csv"), DataFrame)
                    # La función en Resultados.jl maneja la normalización de nombres de columnas de tecnologias_csv_df
                    Resultados.guardar_produccion_diaria_ultimo_dia(df_prod, tecnologias_csv_df, esc_obj)

                    # --- Exportar precios por vector ---
                    try
                        if haskey(res, :precios) && !isempty(res.precios)
                            println("   Exportando precios por vector para $esc_nombre...")
                            Resultados.exportar_precios_vectores(res.precios, esc_obj)
                            
                            df_tec_map = CSV.read("Datos/tecnologias.csv", DataFrame)
                            select!(df_tec_map, [:tecnologia, :vector])
                            # res.opt_resultados[:produccion_mensual] y res.opt_resultados[:deficit_mensual]
                            Resultados.exportar_produccion_mensual_por_vector(
                                res.opt_resultados[:produccion_mensual], df_tec_map, esc_obj
                            )

                            # Exportación **anual** añadida
                            Resultados.exportar_produccion_anual_por_vector(
                                res.opt_resultados[:produccion_mensual],
                                df_tec_map, esc_obj
                            )

                            Resultados.exportar_deficit_mensual(
                                res.opt_resultados[:deficit_mensual], esc_obj
                            )

                        else
                            @warn "   No hay datos de precios disponibles para exportar para $esc_nombre."
                        end
                    catch e
                        @error "Error exportando precios para $esc_nombre: $e"
                        showerror(stdout, e, catch_backtrace()); println()
                        global errores_post += 1
                    end
                    # --- Fin exportación de precios ---
                else
                    @warn "   df_prod vacío para $esc_nombre, omitiendo guardar_produccion_diaria_ultimo_dia."
                end
            catch e
                @error "Error en guardar_produccion_diaria_ultimo_dia para $esc_nombre: $e"
                showerror(stdout, e, catch_backtrace()); println()
                global errores_post += 1
            end
            # --- Fin de nueva llamada ---

        else
             @warn "Índice de resultado $i fuera de rango. No se puede obtener Escenario original para reportes adicionales."
             global errores_post += 1
        end

    elseif opt_status in estados_con_solucion && !isempty(df_cap) # Hay solución subóptima
        @warn "Optimización NO ÓPTIMA para escenario '$esc_nombre', pero con resultados parciales (Estado: $opt_status)."
        opt_obj_val = get(res, :opt_objective, NaN)
        println("   (Valor Objetivo Subóptimo: $(round(opt_obj_val, digits=0)))")
        # Opcional: Guardar/graficar resultados subóptimos
        # ... (similar a arriba, pero quizás con un sufijo "_suboptimo")
        # AÚN ASÍ, intentamos generar los reportes financieros si hay datos
        # (La producción del último día podría no ser relevante si la optimización no fue completa)
        if i <= length(escenarios_para_ejecutar)
            esc_obj_subopt = escenarios_para_ejecutar[i][:escenario_obj] # Renombrar para evitar conflicto de scope si es necesario
            # --- Nueva llamada a exportar_finanzas_agentes_anual (para casos subóptimos también) ---
            try
                finanzas_df_to_export = nothing
                # 1️⃣ Intentamos siempre primero el diccionario anual completo
                if haskey(res, :finanzas_anuales)
                    println("   Usando :finanzas_anuales para reporte financiero de $esc_nombre (subóptimo).")
                    filas = Vector{NamedTuple}()
                    for (_, v) in res.finanzas_anuales
                        append!(filas, v)
                    end
                    finanzas_df_to_export = DataFrame(filas)
                # 2️⃣ Sólo si no existe, recurrimos al antiguo snapshot
                elseif haskey(res, :finanzas_salida) && isa(res.finanzas_salida, DataFrame)
                    finanzas_df_to_export = res.finanzas_salida
                    println("   Usando :finanzas_salida para reporte financiero de $esc_nombre (subóptimo).")
                # 3️⃣ Último recurso: construir la foto del último año desde los agentes
                elseif haskey(res, :agentes)
                    println("   Generando snapshot financiero desde :agentes para $esc_nombre (subóptimo).")
                    ultimo_anio = maximum(keys(res.finanzas_anuales))   # seguro porque existe el dicc.
                    finanzas_df_to_export =
                        DataFrame(SimulacionSecuencial.construir_snapshot_finanzas(res.agentes,
                                                                                    ultimo_anio))
                end

                if finanzas_df_to_export !== nothing && !isempty(finanzas_df_to_export)
                    Resultados.exportar_finanzas_agentes_anual(finanzas_df_to_export, esc_obj_subopt)
                else
                    @warn "   No hay datos financieros disponibles para exportar para $esc_nombre (subóptimo)."
                end
            catch e
                @error "Error en exportar_finanzas_agentes_anual para $esc_nombre (subóptimo): $e"
                showerror(stdout, e, catch_backtrace()); println()
                global errores_post += 1
            end
            # --- Fin de nueva llamada ---
        else
            @warn "Índice de resultado $i fuera de rango (subóptimo). No se puede obtener Escenario original para reportes financieros."
            global errores_post += 1
        end

    else # Estado no exitoso y sin solución o error desconocido
        @warn "Optimización NO EXITOSA para escenario '$esc_nombre' (Estado: $opt_status) y/o sin resultados."
        global errores_post += 1
    end

    # --- Llamada general a exportar_finanzas_agentes_anual para escenarios OK ---
    # Esta sección se moverá/duplicará si queremos que se ejecute también para escenarios óptimos.
    # Por ahora, la he colocado también en el bloque de "subóptimos" si tienen datos.
    # Para mantenerla separada y clara para escenarios óptimos/aceptables:
    if (opt_status in estados_exitosos || es_cournot_ok) && i <= length(escenarios_para_ejecutar)
        esc_obj_opt = escenarios_para_ejecutar[i][:escenario_obj]
        try
            finanzas_df_to_export = nothing

            # 1️⃣ Diccionario anual completo (prioritario)
            if haskey(res, :finanzas_anuales)
                println("   Usando :finanzas_anuales para reporte financiero de $esc_nombre.")
                filas = Vector{NamedTuple}()
                for (_, v) in res.finanzas_anuales
                    append!(filas, v)
                end
                finanzas_df_to_export = DataFrame(filas)

            # 2️⃣ Compatibilidad versiones antiguas
            elseif haskey(res, :finanzas_salida) && isa(res.finanzas_salida, DataFrame)
                finanzas_df_to_export = res.finanzas_salida
                println("   Usando :finanzas_salida para reporte financiero de $esc_nombre.")
                
            elseif haskey(res, :agentes) # Fallback usando agentes finales
                println("   Generando datos financieros desde :agentes para reporte de $esc_nombre.")
                # Obtener último año de la simulación
                anios_simulacion = sort(unique(Int.(datos_bundle_main.demanda.anio)))
                ultimo_anio = maximum(anios_simulacion)
                
                # Crear DataFrame desde agentes finales (sin función auxiliar por ahora)
                filas_agentes = [
                    (
                        agent_id = ag.id,
                        anio = ultimo_anio,
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
                    for ag in values(res.agentes)
                ]
                finanzas_df_to_export = DataFrame(filas_agentes)
            end
            
            # Exportar si hay datos
            if finanzas_df_to_export !== nothing && !isempty(finanzas_df_to_export)
                Resultados.exportar_finanzas_agentes_anual(finanzas_df_to_export, esc_obj_opt)
            else
                @warn "   No hay datos financieros disponibles para exportar para $esc_nombre."
            end
            
        catch e
            @error "Error en exportar_finanzas_agentes_anual para $esc_nombre: $e"
            showerror(stdout, e, catch_backtrace()); println()
            global errores_post += 1
        end
    end
    # --- Fin de llamada general ---


    # Déficit mensual
    rows_def = [(a, m, vec, mwh)
                for (a,dicMes)   in SimulacionSecuencial.deficits_mensual
                for (m,dicVec)   in dicMes
                for (vec,mwh)    in dicVec]
    df_def  = DataFrame(anio = getindex.(rows_def,1),
                        mes  = getindex.(rows_def,2),
                        vector = getindex.(rows_def,3),
                        deficit_mwh = getindex.(rows_def,4))
    CSV.write("deficit_mensual.csv", df_def)


    if !isempty(df_def)
        @df df_def heatmap(:mes, :anio, :deficit_mwh;
                           group = :vector, cbar = true,
                           title = "Déficit (MWh) – heatmap")
    else
        @warn "df_def está vacío; omitiendo heatmap de déficit"
    end

    # Procesar/guardar resultados de otros módulos si es necesario...

end

println("\n" * "-"^40)
println("Resumen de Procesamiento de Resultados:")
total_resultados = length(resultados_completos)
println("  Resultados recibidos: $total_resultados")
println("  Errores durante ejecución en workers (ya logueados): $(
     count(r -> (r isa NamedTuple) && !isnothing(r.error_info),
           resultados_completos)
 )")
println("  Errores durante post-procesamiento (guardar/graficar): $errores_post")
println("✅ Escenarios con optimización ÓPTIMA/ACEPTABLE: $resultados_ok")
println("-"^40)
println("🎉 Ejecución finalizada. Revisa la carpeta '$path_salidas' para resultados detallados y gráficas.")
