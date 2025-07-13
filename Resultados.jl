# =========================================
# Resultados.jl - Exportación y visualización de resultados
# =========================================

module Resultados

using Statistics

using CSV, DataFrames, Plots
gr()    
using StatsPlots # Asegurar StatsPlots para groupedbar
using ..Escenarios: Escenario # Para usar el tipo Escenario
using ..Tecnologias: Tecnologia # Para usar el tipo Tecnologia
# Nota: No necesitamos acceder a Main, pasamos lo necesario como argumentos.

export guardar_resultados,
        graficar_resultados,
        graficar_snapshots,
        graficar_traza,
        append_snapshot!,
        exportar_escenario,
        guardar_produccion_diaria_ultimo_dia,
        exportar_finanzas_agentes_anual,
        exportar_precios_vectores,
        exportar_produccion_anual_por_vector,
        grafico_impacto_negativo_balanza,
        grafico_contaminacion_total_anual

const df_snapshots = DataFrame(
    escenario=String[], anio=Int[], mes=Int[],
    tecnologia=String[],        # <-- nueva
    vector=String[],     # <-- nueva
    capacidad=Float64[], generacion=Float64[],
    emisiones=Float64[], subvencion=Float64[], coste_co2=Float64[]
)



# ═════════════════ Helpers comunes ════════════════════════════════════
# Carga el CSV de tecnologías una sola vez y devuelve sólo :tecnologia ‖ :vector
function _cargar_mapeo_tecnologias()
    tec_map = CSV.read(joinpath("Datos","tecnologias.csv"), DataFrame)
    names!(tec_map, Symbol.(names(tec_map)))          # todo a Symbol
    select!(tec_map, [:tecnologia, :vector])
    return tec_map
end


"""
exportar_produccion_anual_por_vector(df::DataFrame, df_tec::DataFrame, escena)

Genera CSV y gráficas absolutas y relativas de **producción anual** por vector.
Resultados en `Salidas/Distribucion_produccion/<escenario>/`.
"""
function exportar_produccion_anual_por_vector(df::DataFrame,
                                              df_tec::DataFrame,
                                              escena)
    println(" Exportando producción anual para $(escena.nombre)…")
    nombre_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_dp = joinpath("Salidas", "Distribucion_produccion", nombre_safe)
    mkpath(path_dp)

    # df viene en términos mensuales (anio, mes, tecnologia, valor)
    df_vec = leftjoin(df, df_tec, on = :tecnologia)

    # ── CSV total anual por vector ───────────────────────────────────
    df_anual_vec = combine(groupby(df_vec, [:anio, :vector]),
                           :valor => sum => :produccion_mwh)
    CSV.write(joinpath(path_dp,
                       "produccion_anual_por_vector_$(nombre_safe).csv"),
              df_anual_vec)
    # ── Serie temporal absoluta (incluye TOTAL) ─────────────────────
    # 1. Cálculo de la producción TOTAL anual
    df_total = combine(groupby(df_anual_vec, :anio),
                       :produccion_mwh => sum => :produccion_mwh)
    df_total.vector = fill("total", nrow(df_total))        # etiqueta
              
    # 2. Unimos y ordenamos cronológicamente
    df_plot = vcat(df_anual_vec, df_total)
    sort!(df_plot, [:anio])
              
    # 3. Trazado en líneas
    plt_abs = @df df_plot plot(
              :anio, :produccion_mwh;
              group   = :vector,
              lw      = 2,
              marker  = :auto,
              xlabel  = "Año",
              ylabel  = "MWh",
              title   = "Producción anual por vector – $(escena.nombre)")
    savefig(joinpath(path_dp,
                     "produccion_anual_por_vector_$(nombre_safe).png"))
    # ── Relativa (share por tecnología dentro de cada vector) ───────
    df_det = combine(groupby(df_vec, [:anio, :vector, :tecnologia]),
                     :valor => sum => :mwh)
    df_rel = transform(groupby(df_det, [:anio, :vector]),
                       :mwh => (x -> x ./ sum(x)) => :share)
    for v in unique(df_rel.vector)
        df_v = filter(r -> r.vector == v, df_rel)
        plt_rel = @df df_v groupedbar(
            :anio, :share, group = :tecnologia,
            xlabel = "Año", ylabel = "%",
            title  = "Participación relativa $v – $(escena.nombre)")
        savefig(joinpath(path_dp,
                         "participacion_relativa_$(v)_$(nombre_safe).png"))
    end
end

"""
Guardar resultados de la optimización (capacidad, producción) en ficheros CSV.
Argumentos:
- df_capacidad: DataFrame con resultados de capacidad (:tecnologia, :anio, :cap_total, :cap_nueva).
- df_produccion: DataFrame con resultados de producción horaria (:tecnologia, :anio, :estacion, :hora, :valor).
- esc: Objeto Escenario para nombrar los archivos.
"""
function guardar_resultados(df_capacidad::DataFrame, df_produccion::DataFrame, escena) #Antes era esc::Main.Escenarios.Escenario
    
    
    
    # Verificación más robusta
    if :tecnologia in Symbol.(names(df_capacidad))
        rename!(df_capacidad, :tecnologia => :tecnologia)
    end
    
    if :valor in Symbol.(names(df_capacidad))
        rename!(df_capacidad, :valor => :cap_total)
    end

    if :tecnologia in Symbol.(names(df_produccion))
        rename!(df_produccion, :tecnologia => :tecnologia)
    end
        
    # (df_produccion ya tiene :valor y :anio)

    
    
    # Verificar si los DataFrames contienen datos
    if isempty(df_capacidad) && isempty(df_produccion)
        @warn "No hay datos de capacidad ni producción para guardar en escenario $(escena.nombre)."
        return
    end

    println(" Guardando resultados de optimización para $(escena.nombre)...")
    path_salidas = "Salidas"
    isdir(path_salidas) || mkdir(path_salidas) # Crear carpeta si no existe
    # Limpiar archivos PNG anteriores
    println("🧹 Limpiando archivos PNG anteriores...")
    try
        for archivo in readdir(path_salidas; join=true)
            if endswith(lowercase(archivo), ".png")
                rm(archivo, force=true)
                println("   Eliminado: $(basename(archivo))")
            end
        end
        # También limpiar subcarpetas comunes
        subcarpetas = ["DistribucionProduccion_DiaFinal", "Distribucion_produccion", "FinanzasAgentes"]
        for subcarpeta in subcarpetas
            ruta_sub = joinpath(path_salidas, subcarpeta)
            if isdir(ruta_sub)
                for archivo in readdir(ruta_sub; join=true)
                    if isfile(archivo) && endswith(lowercase(archivo), ".png")
                        rm(archivo, force=true)
                    end
                end
            end
        end
        println("✅ Limpieza de PNG completada.")
    catch e
        @warn "Error limpiando archivos PNG: $e"
    end
    # Carpeta por escenario para distribución de producción
    nombre_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_dp = joinpath(path_salidas, "Distribucion_produccion", nombre_safe)
    isdir(path_dp) || mkpath(path_dp)

    # Asegurar columna :agente en los DataFrames de resultados
    # ── Normalizar columna de identificación de agente ───────────
    if :agent_id in names(df_capacidad) && !(:agente in names(df_capacidad))
        rename!(df_capacidad, :agent_id => :agente)
    elseif !(:agente in names(df_capacidad))
        df_capacidad[!, :agente] = fill(missing, nrow(df_capacidad))
    end
    
    if :agent_id in names(df_produccion) && !(:agente in names(df_produccion))
        rename!(df_produccion, :agent_id => :agente)
    elseif !(:agente in names(df_produccion))
        df_produccion[!, :agente] = fill(missing, nrow(df_produccion))
    end

    try
        # Guardar capacidad si no está vacío
        if !isempty(df_capacidad)
            # Validar columnas mínimas antes de guardar
            if all(c -> c in Symbol.(names(df_capacidad)), [:agente, :tecnologia, :anio, :cap_total])
                CSV.write(joinpath(path_salidas, "capacidades_$(escena.nombre).csv"), df_capacidad)
            else
                @warn "Faltan columnas en df_capacidad para $(escena.nombre). No se guardó."
            end
        end

        # Guardar producción (resumen anual) si no está vacío
        if !isempty(df_produccion)
            # Validar columnas mínimas (ya renombradas)
             if all(c -> c in Symbol.(names(df_produccion)), [:agente, :tecnologia, :anio, :valor])
                try
                    # Calcular producción anual
                    df_prod_anual = combine(groupby(df_produccion, [:tecnologia, :anio]), :valor => sum => :produccion_anual_mwh)
                    # Asegurar que el resultado no está vacío
                    if !isempty(df_prod_anual)
                        CSV.write(joinpath(path_salidas, "produccion_anual_$(escena.nombre).csv"), df_prod_anual)
                        
                        # --- CSV y gráficas por vector ---
                        try
                            tec_map = CSV.read(joinpath("Datos","tecnologias.csv"), DataFrame)
                            rename!(tec_map, "tecnologia" => "tecnologia")
                            select!(tec_map, [:tecnologia, :vector])
                            df_prod_vec = leftjoin(df_prod_anual, tec_map, on = :tecnologia)
                            # CSV: producción por vector y tecnología
                            CSV.write(joinpath(path_dp, "produccion_anual_por_vector_$(escena.nombre).csv"), df_prod_vec)
                            # CSV: producción total por vector
                            df_vec_anual = combine(groupby(df_prod_vec, [:vector, :anio]), :produccion_anual_mwh => sum => :produccion_anual_mwh)
                            CSV.write(joinpath(path_dp, "produccion_anual_vector_total_$(escena.nombre).csv"), df_vec_anual)
                            # CSV: participación tecnología dentro de vector
                            df_dist = transform(groupby(df_prod_vec, [:vector, :anio]), :produccion_anual_mwh => (x->round.(x ./ sum(x), digits=4)) => :share)
                            CSV.write(joinpath(path_dp, "participacion_tecnologia_por_vector_$(escena.nombre).csv"), df_dist)
                            # Gráficas por vector
                            for v in unique(df_vec_anual.vector)
                                # Datos totales y por tecnología para vector v
                                df_v = filter(row->row.vector==v, df_vec_anual)
                                df_v_prod = filter(row->row.vector==v, df_prod_vec)
                                
                                # AÑADIR: Asegurar que los datos estén ordenados y sin duplicados
                                df_v = sort(df_v, :anio)
                                df_v = combine(groupby(df_v, :anio), :produccion_anual_mwh => sum => :produccion_anual_mwh)
                                
                                # Gráfico absoluto: total y productores
                                plt_abs = @df df_v plot(:anio, :produccion_anual_mwh; 
                                    label="Total $(v)", 
                                    title="Producción absoluta $(v) – $(escena.nombre)", 
                                    xlabel="Año", 
                                    ylabel="MWh", 
                                    lw=3,
                                    marker=:circle,  # AÑADIR: marcadores para claridad
                                    markersize=4)
                                
                                # MODIFICAR: Agrupar por tecnología antes de graficar
                                for t in unique(df_v_prod.tecnologia)
                                    df_t = filter(row->row.tecnologia==t, df_v_prod)
                                    # AÑADIR: Agregar por año para evitar duplicados
                                    df_t_agg = combine(groupby(df_t, :anio), 
                                                      :produccion_anual_mwh => sum => :produccion_anual_mwh)
                                    df_t_agg = sort(df_t_agg, :anio)
                                    
                                    # Solo graficar si hay datos significativos
                                    if maximum(df_t_agg.produccion_anual_mwh) > 1e-3
                                        plot!(plt_abs, df_t_agg.anio, df_t_agg.produccion_anual_mwh; 
                                              label=t, 
                                              lw=2, 
                                              linestyle=:dash,
                                              marker=:square,
                                              markersize=3)
                                    end
                                end
                                savefig(plt_abs, joinpath(path_dp, "produccion_absoluta_$(v)_$(escena.nombre).png"))
                                
                                # Relativa (distribución tecnología)
                                df_v_rel = df_prod_vec[df_prod_vec.vector .== v, :]
                                if !isempty(df_v_rel)
                                    # AÑADIR: Agregar antes de calcular shares
                                    df_rel = combine(groupby(df_v_rel, [:anio, :tecnologia]), 
                                                    :produccion_anual_mwh => sum => :produccion_anual_mwh)
                                    df_rel = transform(groupby(df_rel, :anio), 
                                                      :produccion_anual_mwh => (x-> x./sum(x)) => :share)
                                    df_rel = sort(df_rel, [:anio, :tecnologia])
                                    
                                    plt_rel = @df df_rel plot(:anio, :share; 
                                        group=:tecnologia, 
                                        seriestype=:line, 
                                        title="Distribución relativa $(v) – $(escena.nombre)", 
                                        xlabel="Año", 
                                        ylabel="Share", 
                                        lw=2,
                                        legend=:outertopright)
                                    savefig(plt_rel, joinpath(path_dp, "participacion_relativa_$(v)_$(escena.nombre).png"))
                                end
                            end
                        catch e
                            @warn "No se pudo generar distribución por vector: $e"
                        end
                    else
                        @warn "Producción anual calculada está vacía para $(escena.nombre)."
                    end
                catch e
                    @error "No se pudo calcular/guardar producción anual para $(escena.nombre): $e"
                end
            else
                 @warn "Faltan columnas en df_produccion para $(escena.nombre). No se guardó resumen anual."
            end
        end
    catch e
         @error "Error general al guardar archivos CSV para $(escena.nombre): $e"
         showerror(stdout, e, catch_backtrace()); println()
    end
end



"""
Graficar y exportar exclusivamente los snapshots (traza mensual).
"""
function graficar_snapshots(escena::String)
    # 0) Traer el global
    df_snap = df_snapshots
    # 1) Volcar CSV de snapshots
    CSV.write("Salidas/snapshots_$(escena).csv", df_snap)

    # 2) Producción total por vector (anual)
    df_vec = combine(groupby(df_snap, [:vector, :anio]),
                     :generacion => sum => :produccion_anual_mwh)
    CSV.write("Salidas/produccion_anual_vector_total_$(escena).csv", df_vec)

    # 3) Producción por vector-tecnología (si existe :tecnologia)
    if :tecnologia in names(df_snap)
        df_det = combine(groupby(df_snap, [:vector, :anio, :tecnologia]),
                         :generacion => sum => :produccion_mwh)
        CSV.write("Salidas/produccion_anual_por_vector_$(escena).csv", df_det)

        # 4) Gráficas por vector
        for v in unique(df_vec.vector)
            # absoluto
            df_tot = filter(r->r.vector==v, df_vec)
            plt_abs = @df df_tot plot(
                :anio, :produccion_anual_mwh;
                xlabel="Año", ylabel="mwh",
                title="Producción absoluta $v – $escena",
                legend=false
            )
            savefig(plt_abs, "Salidas/produccion_absoluta_$(v)_$(escena).png")

            # relativo
            df_vdet = filter(r->r.vector==v, df_det)
            df_rel  = transform(groupby(df_vdet, :anio),
                                :produccion_mwh => (x-> x ./ sum(x)) => :share)
            plt_rel = @df sort(df_rel, [:anio, :tecnologia]) plot(
                :anio, :share;
                group=:tecnologia,
                xlabel="Año", ylabel="Share",
                title="Participación relativa $v – $escena"
            )
            savefig(plt_rel, "Salidas/participacion_relativa_$(v)_$(escena).png")
        end
    end
end








"""
Generar gráficos básicos de resultados de optimización (Capacidad y Producción Anual).
Argumentos:
- resultado_escenario: La NamedTuple devuelta por ejecutar_escenario_worker (contiene .opt_resultados).
- resultado_sim_dummy: Placeholder (no usado actualmente).
- esc: Objeto Escenario original.
- tecnologias: Dict{String, Tecnologia} cargado en el proceso principal.
"""
function graficar_resultados(resultado_escenario::NamedTuple, escena, tecnologias::AbstractDict) # Antes era  esc::Main.Escenarios.Escenario, tecnologias::Dict{String, Tecnologia}
    
    Plots.gr()
    
    println(" Generando gráficos para $(escena.nombre)...")
    path_salidas = "Salidas"
    isdir(path_salidas) || mkdir(path_salidas)

    # Extraer el diccionario de resultados de optimización
    resultados_opt = get(resultado_escenario, :opt_resultados, Dict())
    if isempty(resultados_opt)
         @warn "El campo :opt_resultados está vacío o no existe en el resultado para $(escena.nombre). No se pueden generar gráficos."
         return
    end

    # Extraer DataFrames específicos del diccionario, usando get con default DataFrame()
    df_cap = get(resultados_opt, :capacidad, DataFrame())
    df_prod = get(resultados_opt, :produccion, DataFrame())
    # Podríamos extraer :almacenamiento, :vertido, :reduccion si quisiéramos graficarlos

    # Verificar si tenemos datos para graficar
    if isempty(df_cap) && isempty(df_prod)
        @warn "No hay datos de capacidad ni producción en :opt_resultados para graficar en escenario $(escena.nombre)."
        return
    end

    try
        # --- Gráfico 1: Capacidad total instalada por año y tecnología ---
        if !isempty(df_cap) && all(c -> c in Symbol.(names(df_cap)), [:tecnologia, :anio, :cap_total])
            # Filtrar tecnologías con capacidad > umbral pequeño para evitar ruido
            df_cap_filtrado = filter(row -> row.cap_total > 1e-3, df_cap)

            if !isempty(df_cap_filtrado)
                tecs_con_capacidad = unique(df_cap_filtrado.tecnologia)
                println("   Graficando capacidad para tecnologías: $(join(tecs_con_capacidad, ", "))")

                df_cap_plot = sort(df_cap_filtrado, [:anio, :tecnologia]) # Ordenar para groupedbar

                plt1 = groupedbar(df_cap_plot.anio, df_cap_plot.cap_total, group=df_cap_plot.tecnologia,
                            title="Capacidad Instalada - $(escena.nombre)",
                            xlabel="Año", ylabel="Capacidad Total (mw)",
                            legend=:outertopright, bar_width=0.7, lw=0.1, # line width > 0 para borde
                            framestyle=:box,
                            size=(900, 500)) # Tamaño ajustado
                 savefig(plt1, joinpath(path_salidas, "grafico_capacidad_$(escena.nombre).png"))
                 println("   Gráfico de capacidad guardado.")
            else
                 @warn "No hay capacidad instalada significativa (> 0.001 mw) para graficar en $(escena.nombre)."
            end
        elseif !isempty(df_cap) # Si no está vacío pero faltan columnas
             @warn "DataFrame de capacidad existe pero faltan columnas (:tecnologia, :anio, :cap_total) para $(escena.nombre). Columnas presentes: $(names(df_cap))"
        # else: df_cap está vacío, ya se advirtió antes
        end

        # --- Gráfico 2: Producción anual agregada por tecnología (Stacked Bar) ---
        if !isempty(df_prod) && all(c -> c in Symbol.(names(df_prod)), [:tecnologia, :anio, :valor])
             # Calcular producción anual
             df_prod_anual = combine(groupby(df_prod, [:tecnologia, :anio]), :valor => sum => :prod_anual_mwh)

             # Filtrar tecnologías con producción > umbral
             df_prod_filtrado = filter(row -> row.prod_anual_mwh > 1e-3, df_prod_anual)


            # --- LCOE ---------------------------------------------------------
            df_lcoe = get(resultados_opt, :lcoe, DataFrame())
            if !isempty(df_lcoe)
                plt_lcoe = @df sort(df_lcoe, [:anio, :tecnologia]) plot(
                    :anio, :lcoe_eur_mwh;           #  ← plot en lugar de line
                    group      = :tecnologia,
                    seriestype = :line,             #  ← tipo de serie
                    xlabel     = "Año",
                    ylabel     = "LCOE (€/mwh)",
                    title      = "Coste nivelado – $(escena.nombre)",
                    legend     = :outertopright,
                    framestyle = :box,
                    lw         = 2,
                )
                savefig(plt_lcoe, joinpath(path_salidas,
                        "grafico_LCOE_$(escena.nombre).png"))
                println("   Gráfico LCOE guardado.")
            end

            # --- Gráfico 4: Participación ------------------------------------------
            df_part = get(resultados_opt, :participacion, DataFrame())
            if !isempty(df_part)
                plt_part = @df sort(df_part, [:anio, :tecnologia]) plot(
                :anio, :share;
                group      = :tecnologia,
                seriestype = :line,
                fillrange  = 0,
                xlabel     = "Año",
                ylabel     = "Participación",
                title      = "Adopción por tecnología – $(escena.nombre)",
                legend     = :outerright,
                framestyle = :box,
            )


                savefig(plt_part, joinpath(path_salidas,
                        "grafico_participacion_$(escena.nombre).png"))
                println("   Gráfico de participación guardado.")
            end
            #############

            if !isempty(df_prod_filtrado)
                 tecs_con_produccion = unique(df_prod_filtrado.tecnologia)
                 println("   Graficando producción anual para tecnologías: $(join(tecs_con_produccion, ", "))")

                 df_prod_plot = sort(df_prod_filtrado, [:anio, :tecnologia])

                 plt2 = groupedbar(df_prod_plot.anio, df_prod_plot.prod_anual_mwh, group=df_prod_plot.tecnologia,
                             bar_position = :stack, # Apilar barras
                             title="Producción Anual Agregada - $(escena.nombre)",
                             xlabel="Año", ylabel="Producción Total (mwh)",
                             legend=:outertopright, bar_width=0.7, lw=0.1, # line width > 0 para borde
                             framestyle=:box,
                             size=(900, 500))
                 savefig(plt2, joinpath(path_salidas, "grafico_produccion_anual_$(escena.nombre).png"))
                 println("   Gráfico de producción anual guardado.")
            else
                 @warn "No hay producción anual significativa (> 0.001 mwh) para graficar en $(escena.nombre)."
            end
        elseif !isempty(df_prod) # Si no está vacío pero faltan columnas
             @warn "DataFrame de producción existe pero faltan columnas (:tecnologia, :anio, :valor) para $(escena.nombre). Columnas presentes: $(names(df_prod))"
        # else: df_prod está vacío, ya se advirtió antes
        end

        # ▸ Impacto negativo balanza comercial
        #   df_prod no tiene :vector → unimos con el mapping antes de graficar
        try
            df_prod_vec = leftjoin(df_prod, _cargar_mapeo_tecnologias(), on = :tecnologia)
            # Fallback: si alguna tecnología no aparece en el mapping,
            # utiliza la propia tecnología como vector para no romper el flujo
            if :vector ∉ names(df_prod_vec)
                df_prod_vec[!, :vector] = df_prod_vec.tecnologia
            else
                df_prod_vec[!, :vector] = coalesce.(df_prod_vec.vector, df_prod_vec.tecnologia)
            end
            grafico_impacto_negativo_balanza(df_prod_vec, escena.nombre)
        catch e
            @warn "No se pudo generar 'Impacto negativo balanza comercial': $e"
        end

        # ▸ Contaminación total anual por tecnología
        grafico_contaminacion_total_anual(df_snapshots, escena.nombre)

        println("📊 Gráficos básicos generados para $(escena.nombre).")


        if haskey(resultados_opt, :lcoe)
            CSV.write(joinpath(path_salidas, "lcoe_$(escena.nombre).csv"),
                      resultados_opt[:lcoe])
        end
        if haskey(resultados_opt, :participacion)
            CSV.write(joinpath(path_salidas,
                      "participacion_$(escena.nombre).csv"),
                      resultados_opt[:participacion])
        end
        

    catch e
        @error "Error generando gráficos para $(escena.nombre): $e"
         showerror(stdout, e, catch_backtrace()); println()
    end
end

"""
Graficar la evolución de precio medio por ronda de la traza iterativa.
Argumentos:
- path_traza: Ruta al CSV generado por iterador.
"""
function graficar_traza(path_traza::String)
    try
        df = CSV.read(path_traza, DataFrame)
        # Usar una columna de iteración/paso que sí existe
        x_col = if :anio in names(df)
            :anio
        elseif :mes in names(df)
            :mes
        elseif :step in names(df)
            :step
        elseif :iteracion in names(df)
            :iteracion
        else
            names(df)[1]  # Usar la primera columna como fallback
        end
        
        y_col = if :precio_medio in names(df)
            :precio_medio
        elseif :generacion in names(df)
            :generacion
        elseif :capacidad in names(df)
            :capacidad
        else
            names(df)[end]  # Usar la última columna como fallback
        end
        
        plt = plot(df[!, x_col], df[!, y_col];
                   xlabel=string(x_col), ylabel=string(y_col),
                   title="Evolución - $(basename(path_traza))",
                   legend=false, framestyle=:box)
        ruta_png = replace(path_traza, ".csv" => "_convergencia.png")
        savefig(plt, ruta_png)
        println("   Gráfico de convergencia guardado en $ruta_png")
    catch e
        @warn "Error graficando traza iterativa: $e"
    end
end

"""
Añade un registro mensual al DataFrame global de snapshots.
"""
function append_snapshot!(esc::String, a::Int, m::Int, stats::Dict)
    # stats debe incluir también stats[:tecnologia] y stats[:vector]
    push!(df_snapshots, (
        esc, a, m,
        stats[:tecnologia],       # nombre de tecnología
        stats[:vector],    # vector al que pertenece
        stats[:capacidad_total],
        stats[:generacion_total],
        stats[:emisiones_total],
        stats[:subv_total],
        stats[:coste_co2]
    ))
end

"""
Al final de la simulación, vuelca todos los snapshots a CSV.
"""
function exportar_escenario(esc::String)
    graficar_snapshots(esc)
    CSV.write("Salidas/snapshots_$(esc).csv", df_snapshots)
end


# --- Helper Interno para guardar_produccion_diaria_ultimo_dia ---
"""
_calculate_production_last_day_aggregates(df_produccion_horaria, df_tecnologias; anio_target, dia_target)

Calcula los agregados de producción para un día específico.
Retorna una NamedTuple con (detalle_por_tecnologia, agregado_por_vector) o nothing en caso de error crítico.
"""
function _calculate_production_last_day_aggregates(
    df_produccion_horaria::DataFrame, 
    df_tecnologias::DataFrame; 
    anio_target::Int=2050, 
    dia_target::Int=dayofyear(Date(anio_target,12,31))
)
    # Crear copias para no modificar los originales
    df_prod = copy(df_produccion_horaria)
    df_tec_map = copy(df_tecnologias)

    # Validaciones básicas de entrada
    if isempty(df_prod) || isempty(df_tec_map)
        @warn "_calculate_production_last_day_aggregates: Uno de los DataFrames de entrada está vacío."
        return nothing
    end

    # Normalizar nombres de columnas en df_tec_map (mantener nomenclatura del proyecto)
    if :tecnologia in Symbol.(names(df_tec_map))
        rename!(df_tec_map, :tecnologia => :tecnologia)
    end
    if !all(c -> c in Symbol.(names(df_tec_map)), [:tecnologia, :vector])
        @warn "_calculate_production_last_day_aggregates: df_tecnologias debe contener :tecnologia y :vector."
        return nothing
    end
    select!(df_tec_map, [:tecnologia, :vector])

    # Filtrar producción para el año objetivo
    if !(:anio in Symbol.(names(df_prod))) || !(:estacion in Symbol.(names(df_prod)))
        @warn "_calculate_production_last_day_aggregates: df_prod debe contener :anio y :estacion."
        return nothing
    end
    
    df_prod_target_year = filter(row -> row.anio == anio_target, df_prod)
    if isempty(df_prod_target_year)
        @warn "_calculate_production_last_day_aggregates: No hay datos para el año $anio_target."
        return (detalle_por_tecnologia=DataFrame(), agregado_por_vector=DataFrame())
    end

    # ===== CORRECCIÓN PRINCIPAL: Buscar último periodo disponible =====
    df_ultimo_dia = DataFrame()
    
    # ===== Buscar el día exacto solicitado (por defecto 31 dic) =====
    df_ultimo_dia = filter(r -> get(r, :dia, dia_target) == dia_target, df_prod_target_year)
    if isempty(df_ultimo_dia)
        @warn "No se encontró producción para el día $dia_target del año $anio_target"
        return (detalle_por_tecnologia=DataFrame(), agregado_por_vector=DataFrame())
    end
        
    estaciones_disponibles = unique(df_prod_target_year.estacion)
    if !isempty(estaciones_disponibles)
        # Ordenar estaciones y tomar la última (asumiendo orden: invierno, primavera, verano, otono)
        orden_estaciones = ["invierno", "primavera", "verano", "otono"]
        estaciones_ordenadas = filter(est -> est in estaciones_disponibles, orden_estaciones)
            
        if !isempty(estaciones_ordenadas)
            ultima_estacion = last(estaciones_ordenadas)
            println("   Usando última estación disponible: $ultima_estacion")
            df_ultimo_dia = filter(row -> row.estacion == ultima_estacion, df_prod_target_year)
        else
            # Fallback: usar cualquier estación disponible
            ultima_estacion = String(last(sort(estaciones_disponibles)))
            println("   Usando estación fallback: $ultima_estacion")
            df_ultimo_dia = filter(row -> row.estacion == ultima_estacion, df_prod_target_year)
        end
        
        # Si aún no hay datos, intentar con cualquier dato del año
        if isempty(df_ultimo_dia) && !isempty(df_prod_target_year)
            println("   Usando todos los datos disponibles para $anio_target como último periodo")
            df_ultimo_dia = df_prod_target_year
        end
    end
        

    # Verificación final de df_ultimo_dia
    if isempty(df_ultimo_dia)
        @warn "_calculate_production_last_day_aggregates: No hay datos disponibles para el año $anio_target."
        return (detalle_por_tecnologia=DataFrame(), agregado_por_vector=DataFrame())
    end

    println("   Procesando $(nrow(df_ultimo_dia)) filas de datos de producción para agregación...")

    # Unir con df_tecnologias (mantener nomenclatura del proyecto)
    if !(:tecnologia in Symbol.(names(df_ultimo_dia)))
        @warn "_calculate_production_last_day_aggregates: df_ultimo_dia no tiene columna :tecnologia."
        return nothing
    end
    
    df_prod_vector_detalle = leftjoin(df_ultimo_dia, df_tec_map, on = :tecnologia)

    if !(:valor in Symbol.(names(df_prod_vector_detalle)))
        @warn "_calculate_production_last_day_aggregates: df_prod_vector_detalle no tiene columna :valor."
        return nothing
    end
    
    # Verificar si hay filas sin vector y advertir
    if any(ismissing, df_prod_vector_detalle.vector)
        tecnologias_sin_vector = unique(df_prod_vector_detalle[ismissing.(df_prod_vector_detalle.vector), :tecnologia])
        @warn "_calculate_production_last_day_aggregates: Tecnologías sin vector asignado: $tecnologias_sin_vector"
    end
    
    # Filtrar filas donde vector es missing antes de agrupar
    df_prod_vector_detalle_valid_vector = filter(row -> !ismissing(row.vector), df_prod_vector_detalle)
    if isempty(df_prod_vector_detalle_valid_vector)
        @warn "_calculate_production_last_day_aggregates: No hay datos con vector válido para agregar."
        return (detalle_por_tecnologia=df_prod_vector_detalle, agregado_por_vector=DataFrame())
    end

    # Agregar producción por vector
    df_agregado_vector = combine(groupby(df_prod_vector_detalle_valid_vector, :vector), :valor => sum => :produccion_total_mwh)

    println("   ✓ Agregación completada: $(nrow(df_agregado_vector)) vectores procesados")
    return (detalle_por_tecnologia=df_prod_vector_detalle, agregado_por_vector=df_agregado_vector)
end


"""
guardar_produccion_diaria_ultimo_dia(df_produccion_horaria::DataFrame, df_tecnologias::DataFrame, escena)

Guarda la producción horaria agregada por vector y gráficas para el último día del año 2050 de un escenario específico.

Argumentos:
- df_produccion_horaria: DataFrame con producción horaria. Columnas esperadas: :tecnologia, :anio, :estacion, :hora, :valor.
- df_tecnologias: DataFrame con mapeo de tecnologías a vectores. Columnas esperadas: :tecnologia (o :tecnologia), :vector.
- escena: Objeto o estructura con un campo `nombre` para identificar el escenario.
"""
function guardar_produccion_diaria_ultimo_dia(df_produccion_horaria::DataFrame, df_tecnologias::DataFrame, escena)
    println(" Generando reporte de producción diaria para el último día de 2050 para $(escena.nombre)...")

    # --- 1. Calcular agregados usando el helper interno ---
    # Usar valores por defecto anio_target=2050, dia_target=365
    if !(:tecnologia in Symbol.(names(df_tecnologias)))
        if :tecnologia in Symbol.(names(df_tecnologias))
            rename!(df_tecnologias, :tecnologia => :tecnologia)
        else
            @error "df_tecnologias no tiene columna :tecnologia ni :tecnologia"
            return
        end
    end
    resultados_calculados = _calculate_production_last_day_aggregates(df_produccion_horaria, df_tecnologias; anio_target=2050, dia_target=0)
    
    if isnothing(resultados_calculados) || isempty(resultados_calculados.detalle_por_tecnologia)
        @warn "No se pudieron calcular los agregados de producción del último día para $(escena.nombre). No se generarán archivos."
        return
    end

    df_prod_vector_detalle = resultados_calculados.detalle_por_tecnologia
    df_agregado_vector = resultados_calculados.agregado_por_vector
    
    # --- Setup de Paths y Directorios ---
    nombre_escenario_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_salida_base = "Salidas"
    path_directorio_reporte = joinpath(path_salida_base, "DistribucionProduccion_DiaFinal", nombre_escenario_safe)

    # Crear directorio para gráficas
    path_graficas = joinpath(path_directorio_reporte, "Graficas")
    try
        mkpath(path_graficas)
    catch e
        @error "No se pudo crear el directorio de gráficas: $(path_graficas). Error: $e. Se omitirán las gráficas."
        # No retornar, intentar guardar CSV si es posible
    end

    # --- 2. Generar y Guardar Gráficos por Vector (usando df_prod_vector_detalle) ---
    if !isempty(df_prod_vector_detalle) && all(c -> c in names(df_prod_vector_detalle), [:vector, :tecnologia, :valor])
        println(" Generando gráficos de producción del último día por vector y tecnología...")
        # Filtrar df_prod_vector_detalle para eliminar filas con vector missing antes de iterar unique
        df_plot_base = filter(row -> !ismissing(row.vector), df_prod_vector_detalle)

        for v_actual in unique(df_plot_base.vector) 
            df_vector_especifico = filter(row -> row.vector == v_actual, df_plot_base) # Ya filtrado por !ismissing

            if isempty(df_vector_especifico)
                # Esto no debería ocurrir si unique se aplicó a df_plot_base, pero es una salvaguarda
                @warn "No hay datos de producción para el vector $(v_actual) en el último día de 2050 para $(escena.nombre)."
                continue
            end

            # Agregar producción diaria por tecnología para el vector actual
            df_plot_data = combine(groupby(df_vector_especifico, [:tecnologia]), :valor => sum => :produccion_diaria_mwh)
            
            if isempty(df_plot_data)
                @warn "Datos agregados para graficar están vacíos para el vector $(v_actual) en $(escena.nombre)."
                continue
            end
            
            sort!(df_plot_data, :produccion_diaria_mwh, rev=true)

            try
                titulo_plot = "Producción Último Día 2050 - Vector: $(v_actual) - Escenario: $(escena.nombre)"
                plot_obj = @df df_plot_data bar(
                    :tecnologia, :produccion_diaria_mwh,
                    title=titulo_plot, xlabel="Tecnología", ylabel="Producción (mwh)",
                    legend=false, xrotation=45, bottom_margin=5Plots.mm, size=(800,600)
                )
                
                nombre_archivo_plot = "plot_produccion_ultimodia_2050_vector_$(replace(string(v_actual), r"[^a-zA-Z0-9_]" => "_"))_$(nombre_escenario_safe).png"
                path_completo_plot = joinpath(path_graficas, nombre_archivo_plot)
                savefig(plot_obj, path_completo_plot)
                println("   ✓ Gráfico guardado en: $(path_completo_plot)")
            catch e_plot
                @error "No se pudo generar o guardar el gráfico para el vector $(v_actual) en $(escena.nombre). Error: $e_plot"
                showerror(stdout, e_plot, catch_backtrace()); println()
            end
        end
    else
        @warn "df_prod_vector_detalle está vacío o faltan columnas :vector, :tecnologia, o :valor. No se generarán gráficos para $(escena.nombre)."
    end

    # --- 3. Guardar DataFrame agregado por vector en CSV (usando df_agregado_vector) ---
    if isempty(df_agregado_vector)
        @warn "La producción agregada por vector (para CSV) está vacía para el último día de 2050 en $(escena.nombre). No se guardará CSV."
    else
        try
            mkpath(path_directorio_reporte) # Asegura que el directorio para el CSV exista
            nombre_archivo_csv = "produccion_ultimodia_2050_$(nombre_escenario_safe).csv"
            path_completo_csv = joinpath(path_directorio_reporte, nombre_archivo_csv)
            CSV.write(path_completo_csv, df_agregado_vector)
            println("✓ Reporte CSV de producción del último día de 2050 guardado en: $(path_completo_csv)")
        catch e_csv
            @error "No se pudo crear directorio o guardar el archivo CSV en $(path_directorio_reporte). Error: $e_csv"
        end
    end
end


# --- Helper Interno para exportar_finanzas_agentes_anual ---
"""
_calculate_annual_finanzas_by_agent(finanzas_df::DataFrame)

Agrega los datos financieros de los agentes anualmente.
Retorna un DataFrame agregado o nothing en caso de error.
"""
function _calculate_annual_finanzas_by_agent(finanzas_df::DataFrame)
    df_input = copy(finanzas_df) # Trabajar con una copia

    if isempty(df_input)
        @warn "_calculate_annual_finanzas_by_agent: El DataFrame de entrada está vacío."
        return nothing # O un DataFrame vacío
    end

    columnas_presentes = Symbol.(names(df_input))
    columnas_requeridas_agrupacion = [:agent_id, :anio]

    for col in columnas_requeridas_agrupacion
        if !(col in columnas_presentes)
            @warn "_calculate_annual_finanzas_by_agent: Falta la columna requerida '$(col)'."
            return nothing # O un DataFrame vacío
        end
    end

    columnas_financieras_potenciales = [
        :ingresos, :costo_operacion, :costo_inversion, 
        :beneficio_bruto, :impuestos, :beneficio_neto, :flujo_caja
    ]
    columnas_a_sumar = intersect(columnas_financieras_potenciales, columnas_presentes)

    if isempty(columnas_a_sumar)
        @warn "_calculate_annual_finanzas_by_agent: No se encontraron columnas financieras para sumar."
        return nothing # O un DataFrame vacío
    end
    
    # println("   (Helper) Columnas financieras a agregar: $(join(columnas_a_sumar, ", "))") # Comentado para menos verbosidad

    try
        operaciones_suma = [col => sum .=> col for col in columnas_a_sumar]
        df_finanzas_anual = combine(groupby(df_input, columnas_requeridas_agrupacion), operaciones_suma...)
        
        if isempty(df_finanzas_anual)
            @warn "_calculate_annual_finanzas_by_agent: El DataFrame de finanzas anuales resultó vacío después de la agregación."
            # No es un error per se, podría ser que no haya datos para las combinaciones de agente/año
        end
        return df_finanzas_anual
    catch e
        @error "_calculate_annual_finanzas_by_agent: Error durante la agregación anual: $e"
        showerror(stdout, e, catch_backtrace()); println()
        return nothing # O un DataFrame vacío
    end
end


"""
exportar_finanzas_agentes_anual(finanzas_df::DataFrame, escena)

Agrega los datos financieros de los agentes anualmente y los exporta a un archivo CSV.

Argumentos:
- finanzas_df: DataFrame con datos financieros de agentes. Columnas esperadas: :agente_id, :anio, 
               y métricas financieras como :ingresos, :costo_operacion, etc. Puede incluir :mes.
- escena: Objeto Escenario (o similar con campo `nombre`) para nombrar el directorio y archivo de salida.
"""
function exportar_finanzas_agentes_anual(finanzas_df::DataFrame, escena)
    println(" Exportando finanzas anuales de agentes para el escenario $(escena.nombre)...")

    # --- 1. Calcular agregados usando el helper interno ---
    df_finanzas_anual = _calculate_annual_finanzas_by_agent(finanzas_df)

    if isnothing(df_finanzas_anual) || isempty(df_finanzas_anual) # isempty check for cases where helper returns empty DF
        @warn "No se pudieron calcular las finanzas anuales de agentes para $(escena.nombre) o resultaron vacías. No se guardará CSV."
        return
    end
    
    # --- 2. Crear directorio de salida ---
    nombre_escenario_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_salida_base = "Salidas"
    path_directorio_escenario = joinpath(path_salida_base, "FinanzasAgentes", nombre_escenario_safe)

    try
        mkpath(path_directorio_escenario)
    catch e
        @error "No se pudo crear el directorio de salida: $(path_directorio_escenario). Error: $e"
        return # No continuar si no se puede crear el directorio
    end

    # --- 3. Exportar a CSV ---
    nombre_archivo_csv = "finanzas_anuales_agentes_$(nombre_escenario_safe).csv"
    path_completo_csv = joinpath(path_directorio_escenario, nombre_archivo_csv)

    try
        CSV.write(path_completo_csv, df_finanzas_anual)
        println("✓ Finanzas anuales de agentes guardadas en: $(path_completo_csv)")
    catch e
        @error "No se pudo guardar el archivo CSV de finanzas anuales en $(path_completo_csv). Error: $e"
    end
end


"""
exportar_precios_vectores(precios_dict::Dict, escena)

Exporta precios horarios de cada vector a CSV y genera gráficos de precios promedio anuales.
"""
function exportar_precios_vectores(precios_dict::Dict{Tuple{Int,String,String,Int}, Float64}, escena)
    println(" Exportando precios por vector para el escenario $(escena.nombre)...")
    
    # Crear directorio de precios
    nombre_escenario_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_salida_base = "Salidas"
    path_directorio_precios = joinpath(path_salida_base, "Precios", nombre_escenario_safe)
    
    try
        mkpath(path_directorio_precios)
    catch e
        @error "No se pudo crear el directorio de precios: $(path_directorio_precios). Error: $e"
        return
    end
    
    if isempty(precios_dict)
        @warn "Diccionario de precios vacío para $(escena.nombre). No se exportarán precios."
        return
    end
    
    # Convertir diccionario a DataFrame
    filas_precios = []
    for ((anio, vector, estacion, hora), precio) in precios_dict
        push!(filas_precios, (
            anio = anio,
            vector = vector,
            estacion = estacion,
            hora = hora,
            precio_eur_mwh = precio
        ))
    end
    
    if isempty(filas_precios)
        @warn "No hay datos de precios para exportar para $(escena.nombre)."
        return
    end
    
    df_precios = DataFrame(filas_precios)
    
    # 1. Exportar CSV completo de precios horarios
    nombre_archivo_csv = "precios_horarios_$(nombre_escenario_safe).csv"
    path_completo_csv = joinpath(path_directorio_precios, nombre_archivo_csv)
    
    try
        CSV.write(path_completo_csv, df_precios)
        println("✓ Precios horarios guardados en: $(path_completo_csv)")
    catch e
        @error "Error guardando CSV de precios: $e"
        return
    end
    
    # 2. Calcular precios promedio anuales por vector
    try
        df_precios_anuales = combine(groupby(df_precios, [:anio, :vector]), :precio_eur_mwh => mean => :precio_promedio_anual)
        
        # Guardar CSV de precios anuales
        nombre_archivo_anual_csv = "precios_promedio_anuales_$(nombre_escenario_safe).csv"
        path_anual_csv = joinpath(path_directorio_precios, nombre_archivo_anual_csv)
        CSV.write(path_anual_csv, df_precios_anuales)
        println("✓ Precios promedio anuales guardados en: $(path_anual_csv)")
        
        # 3. Generar gráficos por vector
        vectores_disponibles = unique(df_precios_anuales.vector)
        
        for vector_actual in vectores_disponibles
            df_vector = filter(row -> row.vector == vector_actual, df_precios_anuales)
            
            if nrow(df_vector) < 2
                @warn "Pocos datos para graficar vector $(vector_actual). Saltando gráfico."
                continue
            end
            
            try
                titulo_grafico = "Evolución Precio Promedio Anual - Vector: $(vector_actual) - $(escena.nombre)"
                plot_precios = @df sort(df_vector, :anio) plot(
                    :anio, :precio_promedio_anual,
                    title = titulo_grafico,
                    xlabel = "Año",
                    ylabel = "Precio Promedio (€/mwh)",
                    legend = false,
                    linewidth = 2,
                    marker = :circle,
                    markersize = 4,
                    framestyle = :box,
                    size = (800, 500)
                )
                
                nombre_archivo_plot = "precio_anual_$(replace(string(vector_actual), r"[^a-zA-Z0-9_]" => "_"))_$(nombre_escenario_safe).png"
                path_completo_plot = joinpath(path_directorio_precios, nombre_archivo_plot)
                savefig(plot_precios, path_completo_plot)
                println("   ✓ Gráfico de precios guardado: $(path_completo_plot)")
                
            catch e_plot
                @error "Error generando gráfico de precios para vector $(vector_actual): $e_plot"
            end
        end
        
        # 4. Gráfico conjunto de todos los vectores
        try
            plot_conjunto = @df sort(df_precios_anuales, [:anio, :vector]) plot(
                :anio, :precio_promedio_anual,
                group = :vector,
                title = "Evolución Precios Promedio por Vector - $(escena.nombre)",
                xlabel = "Año",
                ylabel = "Precio Promedio (€/mwh)",
                legend = :outertopright,
                linewidth = 2,
                framestyle = :box,
                size = (1000, 600)
            )
            
            nombre_archivo_conjunto = "precios_todos_vectores_$(nombre_escenario_safe).png"
            path_conjunto = joinpath(path_directorio_precios, nombre_archivo_conjunto)
            savefig(plot_conjunto, path_conjunto)
            println("   ✓ Gráfico conjunto de precios guardado: $(path_conjunto)")
            
        catch e_conjunto
            @error "Error generando gráfico conjunto de precios: $e_conjunto"
        end
        
    catch e
        @error "Error calculando precios anuales para $(escena.nombre): $e"
    end
end


"""
exportar_produccion_mensual_por_vector(df::DataFrame, df_tec::DataFrame, escena)

Genera CSV y gráficas absolutas y relativas de producción mensual por vector.
"""
function exportar_produccion_mensual_por_vector(
    df::DataFrame, df_tec::DataFrame, escena
)
    println(" Exportando producción mensual para $(escena.nombre)...")
    nombre_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_mes = joinpath("Salidas", "Distribucion_produccion_mensual", nombre_safe)
    mkpath(path_mes)

    # Unir con vector
    df_vec = leftjoin(df, df_tec, on = :tecnologia)
    CSV.write(joinpath(path_mes, "produccion_mensual_por_vector_$(nombre_safe).csv"), df_vec)

    # Total absoluto mensual (todos los vectores juntos)
    df_tot = combine(groupby(df_vec, [:anio, :mes]),
                     :valor => sum => :produccion_total)
    sort!(df_tot, [:anio, :mes])                          # orden correcto
    plt_tot = @df df_tot plot(
        :anio .+ (:mes .- 1)/12, :produccion_total;
        xlabel="Mes", ylabel="MWh",
        title="Producción total mensual – $(escena.nombre)",
        legend=false
    )
    savefig(joinpath(path_mes, "produccion_mensual_total_$(nombre_safe).png"))

    # Gráficas por vector
    for v in unique(df_vec.vector)
        df_v = filter(row->row.vector==v, df_vec)
        sort!(df_v, [:anio, :mes])
        # Absoluto
        plt_abs = @df df_v plot(
            :anio .+ (:mes .- 1)/12, :valor;
            xlabel="Mes", ylabel="MWh",
            title="Prod. mensual absoluta $v – $(escena.nombre)",
            legend=false
        )
        savefig(joinpath(path_mes, "produccion_mensual_absoluta_$(v)_$(nombre_safe).png"))
        # Relativo
        df_rel = transform(groupby(df_v, [:anio, :mes]),
                           :valor => (x-> x ./ sum(x)) => :share)
        sort!(df_rel, [:anio, :mes])
        plt_rel = @df df_rel plot(
            :anio .+ (:mes .- 1)/12, :share;
            group=:tecnologia,
            xlabel="Mes", ylabel="Share",
            title="Prod. mensual relativa $v – $(escena.nombre)"
        )
        savefig(joinpath(path_mes, "produccion_mensual_relativa_$(v)_$(nombre_safe).png"))
    end
end


"""
exportar_deficit_mensual(df::DataFrame, escena)

Genera CSV y gráfica conjunta de déficit mensual por vector.
"""
function exportar_deficit_mensual(df::DataFrame, escena)
    println(" Exportando déficit mensual para $(escena.nombre)...")
    nombre_safe = replace(escena.nombre, r"[\\/:*?\"<>|]" => "_")
    path_def = joinpath("Salidas", "Deficit_mensual", nombre_safe)
    mkpath(path_def)

    CSV.write(joinpath(path_def,"deficit_mensual_$(nombre_safe).csv"), df)

    # Gráfica conjunta sólo si existe columna :vector no vacía
    if :vector ∈ names(df) && any(!ismissing, df.vector)
        df_plot = filter(r->!ismissing(r.vector), df)
        sort!(df_plot, [:anio,:mes])
        plt_def = @df df_plot plot(
            :anio .+ (:mes .- 1)/12, :valor;
            group     = :vector,
            xlabel    = "Mes", ylabel = "MWh",
            title     = "Déficit mensual – $(escena.nombre)",
            legend    = :outertopright,
            linewidth = 2,
            framestyle = :box)
    else
        plt_def = @df df plot(
            :anio .+ (:mes .- 1)/12, :valor;
            xlabel    = "Mes", ylabel = "MWh",
            title     = "Déficit mensual – $(escena.nombre)",
            legend    = false,
            linewidth = 2,
            framestyle = :box)
    end
    savefig(joinpath(path_def,"deficit_mensual_conjunto_$(nombre_safe).png"))
end

"""
Grafica las importaciones anuales (columnas cuyo vector empieza por
`import_`) y su suma total.  
Genera:
  • *Salidas/impacto_negativo_balanza_comercial_«esc».png*  
  • *Salidas/impacto_negativo_balanza_comercial_«esc».csv*
"""
function grafico_impacto_negativo_balanza(df_prod::DataFrame,
                                          nombre_escenario::AbstractString)

    # ── 1. Determinar columna idónea para filtrar importaciones ─────────
    col_ref = if     :vector     ∈ names(df_prod)  :vector
              elseif :tecnologia ∈ names(df_prod)  :tecnologia
              else
                  @warn "grafico_impacto_negativo_balanza: DataFrame sin :vector ni :tecnologia."
                  return
              end

    df_imp = filter(col_ref => v -> occursin(r"^import_", String(v)), df_prod)
    isempty(df_imp) && return

    # ── 2. Agregados anuales ────────────────────────────────────────────
    df_imp_agg = combine(groupby(df_imp, [:anio, col_ref]), :valor => sum => :mwh)
    rename!(df_imp_agg, col_ref => :grupo)           # nombre neutro
    df_tot      = combine(groupby(df_imp_agg, :anio), :mwh => sum => :mwh_total)

    # ── 3. Salidas (CSV + PNG) ──────────────────────────────────────────
    dir_out = joinpath("Salidas","Impacto_negativo_balanza_comercial")
    mkpath(dir_out)
    CSV.write(joinpath(dir_out,"importaciones_detalle_anual_$(nombre_escenario).csv"), df_imp_agg)
    CSV.write(joinpath(dir_out,"importaciones_total_anual_$(nombre_escenario).csv"),   df_tot)

    gr()
    plt = plot(title   = "Impacto negativo balanza comercial – $(nombre_escenario)",
               xlabel  = "Año",
               ylabel  = "Importaciones (MWh)",
               legend  = :outertopright,
               size    = (900,500))
    for g in unique(df_imp_agg.grupo)
        d = filter(:grupo => ==(g), df_imp_agg)
        plot!(plt, d.anio, d.mwh; label=String(g), lw=2)
    end
    plot!(plt, df_tot.anio, df_tot.mwh_total; label="Total", lw=3, ls=:dash)
    savefig(plt, joinpath(dir_out,"impacto_negativo_balanza_comercial_$(nombre_escenario).png"))
end

"""
Grafica la contaminación anual emitida por cada tecnología y la suma total
—gráfico de barras apiladas.  
Genera:
  • *Salidas/contaminacion_total_anual_«esc».png*  
  • *Salidas/contaminacion_total_anual_«esc».csv*
"""
function grafico_contaminacion_total_anual(df_snap::DataFrame,
                                           nombre_escenario::AbstractString)
    df_esc = filter(:escenario => ==(nombre_escenario), df_snap)
    isempty(df_esc) && return

    # ── 1. Agregados anuales ───────────────────────────────────────────
    df_tec = combine(groupby(df_esc, [:anio, :tecnologia]), :emisiones => sum => :emis)
    df_tot = combine(groupby(df_esc, :anio),               :emisiones => sum => :emis_total)

    # ── 2. Salidas (CSV + PNG) ─────────────────────────────────────────
    dir_out = joinpath("Salidas","Emisiones")
    mkpath(dir_out)
    CSV.write(joinpath(dir_out,"emisiones_anual_tecnologia_$(nombre_escenario).csv"), df_tec)
    CSV.write(joinpath(dir_out,"emisiones_total_anual_$(nombre_escenario).csv"),      df_tot)

    # ── 3. Gráfico apilado + línea total ──────────────────────────────
    gr()
    df_plot = sort(df_tec, [:anio, :tecnologia])
    plt = groupedbar(df_plot.anio, df_plot.emis;
                     group        = df_plot.tecnologia,
                     bar_position = :stack,
                     xlabel       = "Año",
                     ylabel       = "Emisiones (t CO₂eq)",
                     title        = "Emisiones anuales totales – $(nombre_escenario)",
                     legend       = :outertopright,
                     size         = (900,500))
    plot!(plt, df_tot.anio, df_tot.emis_total;
          seriestype=:line, lw=3, ls=:dash, label="Total")
    savefig(plt, joinpath(dir_out,"grafico_emisiones_anual_$(nombre_escenario).png"))
end

end # module Resultados

