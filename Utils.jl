# =========================================
# Utils.jl - Funciones auxiliares del model
# =========================================

module Utils

# Dependencias externas
using DataFrames
include("Escenarios.jl"); using .Escenarios



# NO 'using .Escenarios' here

export crf, construir_politicas_desde_escenario, normalizar_nombre_tecnologia

"""
Capital Recovery Factor (CRF) para amortizar inversiones.
CRF = (r * (1+r)^n) / ((1+r)^n - 1)
Donde:
- r: tasa de descuento (por periodo, usualmente anual)
- n: número de periodos (usualmente vida útil en años)
"""
function crf(r::Float64, n::Int)::Float64
    # Manejar casos borde
    if r < 0 error("Tasa de descuento 'r' no puede ser negativa.") end
    if n <= 0 error("Número de periodos 'n' debe ser positivo.") end
    # Evitar división por cero si r=0 y n es grande, o si (1+r)^n = 1
    if r ≈ 0.0 return 1.0 / n end # Límite cuando r -> 0
    denominator = (1 + r)^n - 1
    if denominator ≈ 0.0
        # Esto puede pasar si (1+r)^n es muy cercano a 1 (e.g., r muy pequeño, n=1)
        # O podría indicar un problema numérico. Devolver NaN o un valor seguro.
        @warn "CRF denominator is close to zero (r=$r, n=$n). Returning NaN."
        return NaN
    end
    return (r * (1 + r)^n) / denominator
end



"""
Genera un DataFrame con precio de CO₂ por año en base al escenario.
Argumentos:
- escenario: Objeto Escenario que contiene precio inicial e incremento.
- anios: Vector de años para los que generar el precio.
Retorna:
- DataFrame con columnas :anio y :precio_co2 (normalizado).
"""
# *** CORRECCIÓN AQUÍ ***
# Referencia el tipo Escenario usando el nombre del módulo directamente,
# asumiendo que Escenarios está disponible en el scope donde Utils es incluido (Main).
function construir_politicas_desde_escenario(escenario::Escenarios.Escenario, anios::Vector{Int})
    if isempty(anios)
        @warn "Lista de años vacía en construir_politicas_desde_escenario. Devolviendo DataFrame vacío."
        return DataFrame(anio=Int[], precio_co2=Float64[])
    end
    # Calcular precios
    precios = [
        escenario.precio_co2 + (a - minimum(anios)) * escenario.incremento_co2_anual #REVISAR
        for a in anios
    ]
    # Crear DataFrame con nombre de columna normalizado
    return DataFrame(anio = anios, precio_co2 = precios)
end



"""
Normaliza el nombre de una tecnología para uso consistente como clave.
"""
function normalizar_nombre_tecnologia(nombre::AbstractString)::String
    return lowercase(strip(String(nombre)))
end



end # module Utils
