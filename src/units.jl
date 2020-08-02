# TODO : Unitful.jl?
module Units

const s = 1              # seconds
const minute = 60s
const ms = 1e-3s         # milliseconds
const Hz = inv(s)        # Frequency
const kHz = 1e3Hz        # kilohertz
const m = 1              # meter
const cm = 0.01m         # centimeter
const cm² = cm^2         # square centimeter
const mL = cm^3
const L = 1e3mL         # liter
const Liter = L
const μL = 1E-6L
const pL = 1E-12L
const mol = 1           # mole
const M = mol / L       # molar (1000 since the SI units is mM)
const mM = 1E-3M        # millimolar
const μM = 1E-6M        # micromolar
const nM = 1E-9M        # nanomolar
const J = 1              # Joule
const kJ = 1e3J         # kilojoule
const A = 1              # ampere
const μA = 1E-6A         # micrpampere
const V = 1              # volt
const mV = 1E-3V         # millivolt
const S = A / V          # Seimens
const mS = 1E-3S         # milliseimens
const F = A * s / V      # Farad (capacitance)
const Farad = F
const μF = 1E-6F         # microfarad
const C = A * s          # Columb

end # module