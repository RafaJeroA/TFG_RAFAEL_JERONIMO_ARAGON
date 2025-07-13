# =========================================
# InversionesPendientes.jl - Estado compartido de inversiones
# =========================================

module InversionesPendientes

export PENDING_INV, limpiar_pendientes!, obtener_pendientes, registrar_pendiente!, aplicar_pendientes!

using Base.Threads: ReentrantLock
# Registro global de inversiones pendientes (mes_listo, tecnología, mw) por agente
const PENDING_INV = Dict{Int, Vector{Tuple{Int,String,Float64}}}()
const PENDING_LOCK = ReentrantLock()

"""
Limpia todas las inversiones pendientes.
"""
function limpiar_pendientes!()
    empty!(PENDING_INV)
end

"""
Obtiene las inversiones pendientes de un agente.
"""
function obtener_pendientes(agent_id::Int)
    return get(PENDING_INV, agent_id, Vector{Tuple{Int,String,Float64}}())
end

"""
Registra una nueva inversión pendiente para un agente.
"""
function registrar_pendiente!(agent_id::Int, mes_listo::Int, tecnologia::String, mw::Float64)
    lock(PENDING_LOCK) do
        if !haskey(PENDING_INV, agent_id)
            PENDING_INV[agent_id] = Vector{Tuple{Int,String,Float64}}()
        end
        push!(PENDING_INV[agent_id], (mes_listo, tecnologia, mw))
    end
end

"""
Aplica las inversiones pendientes que están listas en el mes actual.
Retorna un vector de tuplas (tecnologia, mw) que fueron aplicadas.
"""
function aplicar_pendientes!(agent_id::Int, mes_actual::Int)
    if !haskey(PENDING_INV, agent_id)
        return Tuple{String,Float64}[]
    end
    
    pending = PENDING_INV[agent_id]
    listas = filter(inv -> inv[1] <= mes_actual, pending)
    aplicadas = [(inv[2], inv[3]) for inv in listas]
    
    # Mantener solo las que aún no están listas
    PENDING_INV[agent_id] = filter(inv -> inv[1] > mes_actual, pending)
    
    return aplicadas
end

export obtener_pendientes_global

"""
Devuelve un `Dict{Int,Vector{Tuple{Int,String,Float64}}}` con
 todas las inversiones pendientes de todos los agentes.
"""
obtener_pendientes_global() = deepcopy(PENDING_INV)

end # module