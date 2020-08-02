##################################
### Commonly-used functions
##################################
"x+1 with the same type as x"
p_one(x) = one(x) + x

"1-x with the same type as x"
one_m(x) = one(x) - x

"""
Regular Hill function
"""
hill(x, k = one(x)) = x / (x + k)
hill(x, k, n) = x^n / (x^n + k^n)

"""
Repressive Hill function
"""
hillr(x, k = one(x)) = hill(k , x)
hillr(x, k, n) = hill(k, x, n)

"""
Michaelis-Menten (MM) function
"""
mm(x, k = one(x)) = hill(x, k)

"""
Repressive MM function
"""
mmr(x, k = one(x)) = hillr(x, k)

"""
Logistic sigmoid function.

See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expit.html
"""
expit(x) = mmr(exp(-x))

"""
    exprel(x, em1 = expm1(x))
Returns x / (exp(x)-1) accurate when x is near zero.
See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exprel.html

Note the fraction is the opposite of `scipy.exprel()`
"""
function exprel(x, em1 = expm1(x))
    res = x / em1
    return ifelse(x ≈ zero(x), one(res), res)
end

"""
Signed sqare root
"""
sqrt_s(x) = sign(x) * sqrt(abs(x))

"""
Signed power
"""
pow_s(x, n) = sign(x) * abs(x)^n

"""
Signed Hill function
"""
hill_s(x, k, n) = sign(x * k) * hill(abs(x), abs(k), n)

"""
Signed repressive Hill function
"""
hillr_s(x, k, n) = hill_s(k, x, n)

"Periodic stimulus current with `strength` for a `duty` for every `peroid`"
i_stim(t; period=one(t), duty=zero(t), strength=zero(t)) = ifelse(mod(t, period) < duty, strength, zero(strength))

"""
    nernst(x_o, x_i[, z=1]; v0=VT)

The Nernst potential ``E_N`` : [Wikipedia](https://en.wikipedia.org/wiki/Nernst_equation) across two compartments.
- species `x` with charge (valence) of `z`,
- `x_o` : the concentration in the outer compartment
- `x_i` : the concentration in the inner compartment
- `v0` : thermal voltage, defaults to 26.71 mV (@310K).
"""
nernst(x_o, x_i; v0=VT) = v0 * log(x_o / x_i)
nernst(x_o, x_i, z::Int; v0=VT) = nernst(x_o, x_i; v0 = VT) / z
const eN = nernst

"""
    pmf(ΔΨ, h_i, h_m; v0 = VT)
    pmf(ΔΨ, ΔpH, ; v0 = VT)

The proton motive force (pmf, ``Δp``) from the mitochondrial membrane potential (`ΔΨ`) and both concentrations of proton in the intermembrane space (`h_i`) and the mitochondrial matrix (`h_m`).
Or from `ΔΨ` and pH difference across (`ΔpH`) the inner mitochodnrial membrane (IMM)
"""
pmf(ΔΨ, h_i, h_m; v0 = VT) = ΔΨ + nernst(h_i, h_m; v0 = VT)
pmf(ΔΨ, ΔpH; v0 = VT) = ΔΨ - log(10) * v0 * ΔpH
const Δp = pmf

"""
    nernst_i(v [, z=1]; iv0=iVT)

The inverse result of Nernst function. Also the reciprocal Boltzmann factor by a voltage difference `v` across a membrane.
"""
nernst_i(v; iv0=iVT) = exp(v * iv0)
nernst_i(v, z::Int; iv0=iVT) = nernst_i(z * v; iv0=iVT)
const ieN = nernst_i

"""
    ghk(p, z, x_i, x_o, zvfrt, em1)
    ghk(p, z, vm, x_i, x_o ; v0 = VT)

GHK flux equation
* `p`: the permeability of the membrane for ion `x`
* `z`: the valence of ion `x`
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


# Conversion between pH and concenrtration
"Convert from power (e.g. pH) to concentration (in mM)"
power_to_conc(p) = exp10(-p + 3)
"Convert from concentration (in mM) to power"
conc_to_power(c) = -log10(c) + 3