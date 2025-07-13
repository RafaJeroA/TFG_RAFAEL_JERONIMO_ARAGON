module PoliticaVerde
export fijar_politica!

"""
    fijar_politica!(anio, mes,
                    emisiones_mes, benef_mes,
                    precio_co2_base, subv_max_restante_anual,
                    gen_verde_mes)

Devuelve `Dict{String,Float64}` con la subvención **por MWh** (negativa → ingreso)
para cada tecnología verde.  
Reglas:
  * Se reparte siempre que quede presupuesto y exista generación verde.
  * El saldo anual restante se distribuye a partes iguales entre los
    meses que faltan (incluido el actual).
"""
function fijar_politica!(
        anio::Int, mes::Int,
        _emisiones_mes::Float64,     # (no usado por la regla actual)
        _benef_mes::Dict{Int,Float64},
        _precio_co2_base::Float64,
        subv_max_restante_anual::Float64,
        gen_verde_mes::Dict{String,Float64})

    total_green_mwh = sum(values(gen_verde_mes))

    if total_green_mwh == 0 || subv_max_restante_anual ≤ 0.0
        return Dict{String,Float64}()        # nada que subvencionar
    end

    # Presupuesto para el mes = resto anual / meses pendientes (incluye el actual)
    meses_pendientes = 12 - mes + 1
    presup_mes = subv_max_restante_anual / meses_pendientes

    subv_unit = -presup_mes / total_green_mwh   # negativa ⇒ ingreso agente
    return Dict(k => subv_unit for k in keys(gen_verde_mes))
end


end # module
