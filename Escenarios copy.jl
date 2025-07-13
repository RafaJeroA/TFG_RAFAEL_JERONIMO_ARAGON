# =========================================
# Escenarios.jl - Definición de políticas y shocks
# =========================================

module Escenarios

export Escenario, definir_escenarios, crear_escenario

"""
Estructura de un escenario para simulación y optimización.
Contiene políticas, precios y eventos que afectan la transición.
"""
mutable struct Escenario
    nombre::String
    precio_co2::Dict{Tuple{Int64,Int64}, Float64}
    subv_max_anual::Dict{Int,Float64}   # € máx. para subvención cada año   
    politica_vehiculos_limpios::Bool
    shock_petroleo::Bool
    retiro_programado::Dict{String, Dict{Int, Float64}}             # Campo para gestionar la retirada programada de tecnologías (tec -> año -> MW)
    modo_debug::Bool  
end


"""
Define una lista de escenarios para evaluación comparativa.
"""
function definir_escenarios()
    # Lista de políticas base: (nombre, precio_co2, subv_max, vehículos_limpios, shock_petroleo)
    # vehiculos limpios iba a ser un impuesto extra a em co2 coches y shock petrolero era un limite extricto a importaciones de comb fósiles
    escenarios = [    
        ("Base",           Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), true, false),
        ("BaseSinCierre",           Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), false, false),
        ("BAU",            Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), false, false),
        ("Esfuerzo_x2",    Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), true,  false), 
        ("Esfuerzo_x2SinCierre",    Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), true,  false), 
        ("Esfuerzo_mitad", Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(), true,  false)
    ]
    
    for (nombre, pc, subv, vlimp, shock) in políticas
        # Variante sin cierre nuclear
        push!(escenarios,
             crear_escenario(nombre;
                             precio_co2=pc,
                             subv_max=subv,
                             vehiculos_limpios=vlimp,
                             shock_petroleo=shock,
                             cierre_nuclear=false))
        # Variante con cierre nuclear
        push!(escenarios,
             crear_escenario(nombre;
                             precio_co2=pc,
                             subv_max=subv,
                             vehiculos_limpios=vlimp,
                             shock_petroleo=shock,
                             cierre_nuclear=true))
    end
    return escenarios
end

"""
Constante con los MW que se retiran cada año si hay cierre nuclear PROGRAMACIÓN DEL GOBIERNO
"""
const RETIROS_NUCLEAR = Dict(
    2027 => 1049.40,  # Almaraz I
    2028 => 1044.50,  # Almaraz II
    2030 => 2124.52,  # Ascó I + Cofrentes
    2032 => 1027.21,  # Ascó II
    2025 => 2153.14   # Vandellós II + Trillo
)

"""
Helper para crear un escenario.
- nombre_base: nombre sin sufijo
- precio_co2, subv_max, vehiculos_limpios, shock_petroleo: parámetros habituales
- cierre_nuclear: si true, inyecta RETIROS_NUCLEAR en el dict
Devuelve un Escenario con nombre extendido:
  nombre_base_SinCierre   o   nombre_base_CierreNuclear
"""
function crear_escenario(nombre_base::String;
                          precio_co2=Dict{Tuple{Int64,Int64},Float64}(),
                          subv_max=Dict{Int,Float64}(),
                          vehiculos_limpios::Bool=false,
                          shock_petroleo::Bool=false,
                          cierre_nuclear::Bool=false,
                          modo_debug::Bool=true)
    retiro = cierre_nuclear ? Dict("nuclear" => RETIROS_NUCLEAR) : Dict{String,Dict{Int,Float64}}()
    nombre = cierre_nuclear ? "$(nombre_base)_CierreNuclear" : "$(nombre_base)_SinCierre"
    return Escenario(nombre,
                     precio_co2,
                     subv_max,
                     vehiculos_limpios,
                     shock_petroleo,
                     retiro,
                     modo_debug)
end


end # module Escenarios