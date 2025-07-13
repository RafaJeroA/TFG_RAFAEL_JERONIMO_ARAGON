# =========================================
# Tecnologias.jl - Parámetros y evolución tecnológica
# =========================================

module Tecnologias

using DataFrames, Statistics
export Tecnologia, inicializar_tecnologias!, calcular_LCOE, actualizar_costo_aprendizaje!, actualizar_costo_marginal_aprendizaje!
using ..Utils: normalizar_nombre_tecnologia

using CSV

"""
Estructura que define una tecnología energética en el model.
Incluye parámetros técnicos, económicos y ambientales.
"""
mutable struct Tecnologia
    nombre::String
    tipo::String              # Ej: "generacion", "movilidad", "almacenamiento", "h2"
    vector::String            # Ej: "electricidad", "hidrogeno", "combustibles", SALIDA
    eficiencia::Float64
    vida_util::Int
    emisiones_unitarias::Float64    # tco2 por unidad producida
    tasa_aprendizaje::Float64       # Porcentaje reducción por duplicación de capacidad
    costo_inicial::Float64          # Coste base inicial (€/mw [o mw equivalente])
    costo_om::Float64               # Coste anual de operación (€/mw/anio)
    capacidad_acumulada::Float64    # Capacidad instalada hasta la fecha (mw)
    curva_LCOE::Dict{Int, Float64}  # LCOE anual estimado
    costo_om_variable::Float64 
    vector_consumido::String    # vector de entrada para procesos de conversión (vacío si no aplica)
    ratio_energia::Float64      # mwh producido por mwh consumido en procesos de conversión
    aporta_inercia::Float64
    verde::Bool                    # True si renovable sin CO₂
    tasa_autodescarga::Float64
end

"""
sumar_capacidad_inicial() -> Dict{String,Float64}

Lee los CSV de agentes de generación y de tecnologías, filtra sólo las
tecnologías declaradas en `tecnologias.csv` y suma la capacidad (MW) de
cada una en todos los agentes. Devuelve un Dict con clave = nombre de
tecnología y valor = capacidad total instalada (MW).
"""
function sumar_capacidad_inicial()::Dict{String,Float64}
    # Cargar datos
    agentes_df_inicial    = CSV.read("Datos/agentes_generacion.csv", DataFrame)
    tecnologias_df_inicial = CSV.read("Datos/tecnologias.csv", DataFrame)

    # Conjunto de tecnologías válidas
    tecnos = Set(tecnologias_df_inicial.tecnologia)

    # Filtrar agentes cuya tecnología esté en la lista oficial
    df_filtrado = filter(row -> row.tecnologia in tecnos, agentes_df_inicial)

    # Agrupar por tecnología y sumar capacidad
    cap_inicial_tecs = combine(groupby(df_filtrado, :tecnologia),
                     :capacidad_mw => sum => :capacidad_total_inicial)

    # Convertir a Dict para acceso rápido
    return Dict(row.tecnologia => row.capacidad_total_inicial for row in eachrow(cap_inicial_tecs))
end


const capacidad_inicial_tecnologias = sumar_capacidad_inicial()

"""
Actualiza el **coste marginal variable** (`costo_om_variable`, € /MWh) de
una tecnología aplicando su `tasa_aprendizaje`.  
La reducción se modela igual que en CAPEX: por cada duplicación de la
capacidad acumulada el coste baja un `tasa_aprendizaje` %.

```
factor = (capacidad_acumulada / base)^(-log₂(1 - τ))
nuevo_coste = coste_actual * factor
```

*Parámetros*  
`base` evita divisiones por 0 (valor por defecto = 1 MW).  
La función modifica el campo `costo_om_variable` **in-place** y lo devuelve.
"""
function actualizar_costo_marginal_aprendizaje!(tecnologia::Tecnologia) #VOY A PONER COMO BASE LOS 1000 MW AL CONSIDERAR QUE A PARTIR DE AHÍ EMPIEZA LA ECONOMÍA DE ESCALA Y LO QUE HAY POR ABAJO PUEDE SER UN PROTOTIPO INCLUSO
    base = get(capacidad_inicial_tecnologias, tecnologia.nombre, 100.0)
    cap = max(tecnologia.capacidad_acumulada, base) #BASE = COSTE MARGINAL VARIABLE DEL INICIO
    factor = (cap / base)^(-log2(1 - tecnologia.tasa_aprendizaje))
    tecnologia.costo_om_variable *= factor
    return tecnologia.costo_om_variable
end


"""
Inicializa todas las tecnologías con sus parámetros base desde el DataFrame de datos.
Aplica la curva de aprendizaje si hay capacidad acumulada.
"""
function inicializar_tecnologias!(df::DataFrame)::Dict{String,Tecnologia}

    println("🔍 DEBUG: Columnas del DataFrame tecnologías: $(names(df))")
    println("🔍 DEBUG: Primeras 3 filas de tecnología:")
    if nrow(df) >= 3
        for i in 1:3
            println("  $(i): $(df[i, :tecnologia])")
        end
    end
    
        
        # Resto de la función...
    # Convertir costes de € sobre el DataFrame de tecnologías
    # Convertir costes de forma segura, manejando missing values
    df.costo_inicial = [ismissing(x) ? 0.0 : Float64(x) for x in df.costo_inicial]
    df.costo_om = [ismissing(x) ? 0.0 : Float64(x) for x in df.costo_om]
    if :costo_om_variable in names(df)
        df.costo_om_variable = [ismissing(x) ? 0.0 : Float64(x) for x in df.costo_om_variable]
    end

    println("🔍 DEBUG: Filas después de conversión segura: $(nrow(df))")
    println("🔍 DEBUG: Tecnologías que se van a cargar:")
    for (i, tec) in enumerate(df.tecnologia[1:min(10, nrow(df))])
        println("  $i: $tec")
    end


    tecs = Dict{String,Tecnologia}()
    for (i, row) in enumerate(eachrow(df))
        nombre_original = row.tecnologia
        nombre = normalizar_nombre_tecnologia(nombre_original)
        
        # Debug para ver qué tecnología se está procesando
        if nombre == "coche_gasolina"
            println("🔍 DEBUG: Procesando coche_gasolina en fila $i")
            println("  - costo_inicial: $(row.costo_inicial)")
            println("  - costo_om: $(row.costo_om)")
            println("  - costo_om_variable: $(get(row, :costo_om_variable, "NO_EXISTE"))")
        end

        # Aquí definimos costo_var a partir de la columna, o 0.0 si no existe/missing
        costo_var = coalesce(get(row, :costo_om_variable, 0.0), 0.0)

        try
            tecs[nombre] = Tecnologia(
                nombre,
                row.tipo,
                row.vector,
                row.eficiencia,
                row.vida_util,
                row.emisiones_unitarias,
                coalesce(row.tasa_aprendizaje, 0.0),
                row.costo_inicial,
                row.costo_om,
                0.0,
                Dict{Int, Float64}(),
                costo_var,
                coalesce(get(row, :vector_consumido, ""), ""),
                begin
                    r_tmp = get(row, :ratio_energia, 1.0)
                    # Si es missing, NaN, infinito o ≤0 → valor seguro
                    if ismissing(r_tmp) || !isfinite(r_tmp) || r_tmp ≤ 0
                        #@warn "ratio_energia inválido ($r_tmp) para la tecnología '$(nombre)'; se usa 1.0."   es normal que no haya si no se consume nada
                        r_tmp = 1.0
                    end
                    Float64(r_tmp)
                end,
                coalesce(get(row, :aporta_inercia, 0.0), 0.0),
                (row.emisiones_unitarias == 0.0),
                coalesce(get(row, :tasa_autodescarga, 0.0), 0.0)
            )
            
            if nombre == "coche_gasolina"
                println("✅ DEBUG: coche_gasolina añadida correctamente al diccionario")
            end
        catch e
            @error "Error creando tecnología '$nombre': $e"
            if nombre == "coche_gasolina"
                println("❌ DEBUG: Error creando coche_gasolina: $e")
            end
        end
    end

    println("🔍 DEBUG: Total tecnologías creadas: $(length(tecs))")
    println("🔍 DEBUG: ¿coche_gasolina en diccionario final? $(haskey(tecs, "coche_gasolina"))")
    return tecs
end


"""
Actualiza el coste de inversión de una tecnología aplicando curva de aprendizaje.
Reducción por duplicación: nuevo_coste = coste * (cap_acum / base)^(-log2(1 - tasa))
Se actualiza tecnologia.costo_inicial directamente.
"""
function actualizar_costo_aprendizaje!(tecnologia::Tecnologia, base::Float64=1.0)
    cap = max(tecnologia.capacidad_acumulada, base)  # evitar división por cero
    factor = (cap / base)^(-log2(1 - tecnologia.tasa_aprendizaje))
    tecnologia.costo_inicial *= factor
    return tecnologia.costo_inicial
end



end # module Tecnologias