module ModelShortcodes

using UnPack

export T0, F, R, VT, iVT, KEQ_MGATP, valence, power_to_conc, conc_to_power, Chemical, p_one, one_m, mm, mmr, hill, hillr, expit, nernst, pmf, ghk, sqrt_s, pow_s, exprel, i_stim, hill_s, hillr_s, binding, poly, binding_fraction,
fractions

# TODO: Unitful.jl?
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

using .Units: C, mol, J, kJ, M, mM, μM

# Physical constants
const T0 = 310             # Default temp (37C)
const F = 96485C/mol       # Faraday constant
const R = 8.314J/mol       # Ideal gas constant
const VT = R * T0 / F      # Thermal voltage (@37C)
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)
const ΔG_MGATP = -29.08kJ/mol  # Free energy change for MgATP + H2O <=> MgADP + Pi + H
const KEQ_MGATP = exp(-ΔG_MGATP / (R * T0))   # equlibrium constant of MgATP hydrolysis (37C, pH = 7, [Mg] = 1mM)

# Chemical Tags
abstract type Chemical end
abstract type Cation <: Chemical end
abstract type Anion <: Chemical end
abstract type Phosphate <: Anion end
abstract type Adenylate <: Anion end
abstract type CationV1 <: Cation end
abstract type CationV2 <: Cation end

struct Water <: Chemical end
struct HydrogenPeroxide <: Chemical end
struct Superoxide <: Anion end
struct Hydroxide <: Anion end
struct Sodium <: CationV1 end
struct Potassium <: CationV1 end
struct Proton <: CationV1 end
struct Calcium <: CationV2 end
struct Magnesium <: CationV2 end
struct ATP <: Adenylate end
struct ADP <: Adenylate end
struct AMP <: Adenylate end
struct HATP <: Adenylate end
struct HADP<: Adenylate  end
struct MgATP <: Adenylate end
struct MgADP <: Adenylate end
struct HPO4 <: Phosphate end
struct H2PO4 <: Phosphate end

"Ion valence"
valence(::Type{T}) where {T<:Chemical} = 0
valence(::Type{T}) where {T<:CationV1} = +1
valence(::Type{T}) where {T<:CationV2} = +2
valence(::Type{Superoxide}) = -1
valence(::Type{Hydroxide}) = -1
valence(::Type{ATP}) = -4
valence(::Type{ADP}) = -3
valence(::Type{AMP}) = -2
valence(::Type{Phosphate}) = -3
valence(::Type{HATP}) = valence(ATP) + valence(Proton)
valence(::Type{HADP}) = valence(ADP) + valence(Proton)
valence(::Type{MgATP}) = valence(ATP) + valence(Magnesium)
valence(::Type{MgADP}) = valence(ADP) + valence(Magnesium)
valence(::Type{HPO4}) = valence(Phosphate) + valence(Proton)
valence(::Type{H2PO4}) = valence(Phosphate) + 2 * valence(Proton)
valence(::Type{Water}) = valence(Hydroxide) + valence(Proton)



# Commonly-used functions

"x+1 with the same type as x"
p_one(x) = one(x) + x

"1-x with the same type as x"
one_m(x) = one(x) - x

"Michaelis-Menten (MM) function"
mm(x, k = one(x)) = x / (x + k)

"Repressive MM function"
mmr(x, k = one(x)) = k / (k + x)

"Regular Hill function"
hill(x, k, n) = mm((x / k)^n)

"Repressive Hill function"
hillr(x, k, n) = mmr((x / k)^n)

"Logistic sigmoid function"
expit(x) = mmr(exp(-x))

"Returns x / (exp(x)-1) accurate when x is near zero"
exprel(x, em1 = expm1(x)) = ifelse(x ≈ zero(x), one(x), x / em1)

"Give signs to functions as workarounds for negative values"
_sign(f::Function, x) = sign(x) * f(abs(x))

"Signed sqare root"
sqrt_s(x) = _sign(sqrt, x)

"Signed power"
pow_s(x, n) = _sign(x -> x^n, x)

"Signed Hill function"
hill_s(x, k, n) = _sign(x -> hill(x, k, n), x)

"Signed repressive Hill function"
hillr_s(x, k, n) = _sign(x -> hill_s(x, k, n), x)

"Periodic stimulus current with `strength` for a `duty` for every `peroid`"
i_stim(t; period=one(t), duty=zero(t), strength=zero(t)) = ifelse(mod(t, period) < duty, strength, zero(strength))

"""
    nernst(x_o, x_i [, z=1]; v0=VT)
    nernst(x_o, x_i [, X]; v0=VT)

Get the Nernst potential across two compartments of species `X` with charge (valence) of `z`, where `x_o` is the concentration in the outer compartment, while `x_i` is the inner one.
OPtional parameter `v0` indicates the thermal voltage, defaults to 26.71 mV, corresponding the temprature of 310K.
"""
nernst(x_o, x_i; v0=VT) = v0 * log(x_o / x_i)
nernst(x_o, x_i, z::Int; v0=VT) = nernst(x_o, x_i; v0 = VT) / z
nernst(x_o, x_i, X::Type{T}; v0=VT) where {T<:Chemical} = nernst(x_o, x_i, valence(X); v0=VT)

"""
    pmf(ΔΨ, h_i, h_m; v0 = VT)
    pmf(ΔΨ, ΔpH, ; v0 = VT)

The proton motive force (pmf) from the mitochodnrial membrane potential (`ΔΨ`) and both concentrations of proton in the intermembrane space (`h_i`) and the mitochondrial matrix (`h_m`).
Or from `ΔΨ` and pH difference across (`ΔpH`) the inner mitochodnrial membrane (IMM)
"""
pmf(ΔΨ, h_i, h_m; v0 = VT) = ΔΨ + nernst(h_i, h_m; v0 = VT)
pmf(ΔΨ, ΔpH; v0 = VT) = ΔΨ - log(10) * v0 * ΔpH

"""
    ghk(p, z, x_i, x_o, zvfrt, em1)
    ghk(p, z, vm, x_i, x_o ; v0 = VT)
    ghk(p, X::Chemical, vm, x_i, x_o; v0 = VT)

GHK flux equation
* `p`: the permeability of the membrane for ion `X`
* `z`: the valence of ion `X`
* `v`: the transmembrane potential in volts
* `x_i`: the internal concentration of ion S, measured in mol·m−3 or mM
* `x_o`: the external concentration of ion S, measured in mol·m−3 or mM
* `zvfrt`: precalculated `z * vm / VT`
* `em1`: precalculated `exp(zvfrt) - 1`
* `v0`: the thermal voltage, defaults to 26.71 mV (@ 37C)
"""
ghk(p, z, x_i, x_o, zvfrt, em1) = p * z * F * exprel(zvfrt, em1) * (p_one(em1) * x_i - x_o)

function ghk(p, z, vm, x_i, x_o; v0 = VT)
    zvfrt = z * vm / VT
    em1 = expm1(zvfrt)
    return ghk(p, z, x_i, x_o, zvfrt, em1)
end

ghk(p, X::Type{T}, vm, x_i, x_o; v0 = VT) where {T<:Chemical} = ghk(p, valence(X), vm, x_i, x_o; v0 = VT)

# Conversion between pH and concenrtration
"Convert from power (e.g. pH) to concentration (in mM)"
power_to_conc(p) = exp10(-p + 3)
"Convert from concentration (in mM) to power"
conc_to_power(c) = -log10(c) + 3

# Constants from Guynn 1973 (pH 7.4, 38C, Mg = 1mM, ionic strength = 0.2 M)
binding(X::Type{T}, Y::Type{S}) where {T<:Chemical, S<:Chemical} = Inf  # defaults to infinity
binding(X::Type{T}, Y::Type{S}) where {T<:Chemical, S<:Cation} = binding(Y, X) # Cations could be before or after anions
binding(X::Type{T}, Y::Type{T}) where {T<:Cation} = Inf                # cation to cation: Inf
binding(X::Type{Proton}, Y::Type{Hydroxide}) = 1E-14M^2
binding(X::Type{Proton}, Y::Type{HPO4}) = 1.76E-7M
binding(X::Type{Proton}, Y::Type{ATP}) = 1.08E-7M
binding(X::Type{Proton}, Y::Type{ADP}) = 1.20E-7M
binding(X::Type{Proton}, Y::Type{AMP}) = 3.24E-7M
binding(X::Type{Magnesium}, Y::Type{ATP}) = inv(13900/M)
binding(X::Type{Magnesium}, Y::Type{HATP}) = inv(35.5/M)
binding(X::Type{Magnesium}, Y::Type{ADP}) = inv(1320/M)
binding(X::Type{Magnesium}, Y::Type{HADP}) = inv(32.4/M)
binding(X::Type{Magnesium}, Y::Type{AMP}) = inv(60.1/M)
binding(X::Type{Magnesium}, Y::Type{HPO4}) = inv(90.31/M)

"Binding fraction"
binding_fraction(x, X::Type{T}, Y::Type{S}) where {T<:Chemical, S<:Chemical} = x / binding(X, Y)

"Binding polynomial, sum of binding fractions"
poly(x, X::Type{T}, Y::Type{S}) where {T<:Chemical, S<:Chemical} = p_one(binding_fraction(x, X, Y))

"Adenine nucleotide subcomponents fractions"
function fractions(h, mg, A::Type{T}) where T <: Adenylate
    fH = binding_fraction(h, Proton, A)
    fMg = binding_fraction(mg, Magnesium, A)
    p = p_one(fH + fMg)
    axpn = inv(p)
    return (axpn = axpn, haxp = fH * axpn, mgaxp = fMg * axpn, p = p)
end

"Phosphate to its binding components"
function fractions(h, ::Type{T}) where {T <: Phosphate}
    p = poly(h, Proton, HPO4)
    hpo4 = inv(p)
    return (hpo4 = hpo4, h2po4 = omx(hpo4), p = p)
end

end
