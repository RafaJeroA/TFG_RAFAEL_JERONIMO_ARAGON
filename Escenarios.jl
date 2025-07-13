# =========================================
# Escenarios.jl - Definición de políticas y shocks
# =========================================

module Escenarios

export Escenario, definir_escenarios

# ───────────────── Calendario de cierre nuclear ───────────────────────────
const CALENDARIO_CIERRE_NUCLEAR = Dict(
    2027 => 1_049.40,   # Almaraz I
    2028 => 1_044.50,   # Almaraz II
    2030 => 2_124.52,   # Ascó I + Cofrentes
    2032 => 1_027.21,   # Ascó II
    2035 => 2_153.14    # Vandellós II + Trillo #ES 203V, POR SI LO CAMBIO CON REPLACE SIN QUERER
)


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

    esc_base = Escenario("Base",           Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(),
                         false, false, Dict(), true)
    esc_bau  = Escenario("BAU",            Dict{Tuple{Int,Int},Float64}(),    Dict{Int,Float64}(),
                         false, false, Dict(), false)
    esc_x2   = Escenario("Esfuerzo_x2",    Dict{Tuple{Int,Int},Float64}(),    Dict{Int,Float64}(),
                         false, false, Dict(), false)
    esc_x05  = Escenario("Esfuerzo_mitad", Dict{Tuple{Int,Int},Float64}(),    Dict{Int,Float64}(),
                         false, false, Dict(), false)

    esc_base_sin_cierre = Escenario("BaseSinCierre",           Dict{Tuple{Int64,Int64},Float64}(), Dict{Int,Float64}(),
                         false, false, Dict(), true)
    esc_bau_sin_cierre  = Escenario("BAUSinCierre",            Dict{Tuple{Int,Int},Float64}(),    Dict{Int,Float64}(),
                         false, false, Dict(), false)
    esc_x2_sin_cierre   = Escenario("Esfuerzo_x2SinCierre",    Dict{Tuple{Int,Int},Float64}(),    Dict{Int,Float64}(),
                         false, false, Dict(), false)



    # ─── Inyectar el calendario de cierre en los escenarios solicitados ───
    for esc in [esc_base, esc_bau, esc_x2, esc_x05]
        esc.retiro_programado = Dict("nuclear" => CALENDARIO_CIERRE_NUCLEAR)
    end
     
    return [esc_base, esc_bau, esc_x2, esc_x05, esc_base_sin_cierre, esc_bau_sin_cierre, esc_x2_sin_cierre]
end

end # module Escenarios