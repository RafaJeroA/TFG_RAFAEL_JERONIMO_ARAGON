module FinanzasParametros

export INT_RATE, PLAZO_PRESTAMO, EQ_SHARE, DSCR_MIN, ND_EBITDA_MAX, TAX_RATE, DEP_METHOD

const INT_RATE        = 0.02      # 2 % interés fijo
const PLAZO_PRESTAMO       = 15        # años
const EQ_SHARE        = 0.20      # 20 % equity
const DSCR_MIN        = 1.30      # Cobertura mínima de servicio de deuda
const ND_EBITDA_MAX   = 4.0       # Ratio deuda/EBITDA máximo
const TAX_RATE        = 0.10      # Impuesto sobre EBT
const DEP_METHOD      = :lineal   # Depreciación anual = CAPEX/vida

end
