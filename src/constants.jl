# Physical constants
import .Units: C, mol, J, kJ, M, mM, μM
const T0 = 310             # Default temp (37C)
const F = 96485C/mol       # Faraday constant
const R = 8.314J/mol       # Ideal gas constant
const VT = R * T0 / F      # Thermal voltage (@37C)
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)
const ΔG_MGATP = -29.08kJ/mol  # Free energy change for MgATP + H2O <=> MgADP + Pi + H
const KEQ_MGATP = exp(-ΔG_MGATP / (R * T0)) * M   # equlibrium constant of MgATP hydrolysis (37C, pH = 7, [Mg] = 1mM)
const KW = 1e-14 * M^2     # Ion product of water
