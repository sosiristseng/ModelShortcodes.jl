using .Units: M, mM, μM
# Chemical Tags
abstract type Chemical end
abstract type Cation <: Chemical end
abstract type Anion <: Chemical end
abstract type Phosphate <: Anion end
abstract type Adenylate <: Anion end
abstract type Cation1 <: Cation end
abstract type Cation2 <: Cation end

"Get quantity values chemical object"
qty(x::Chemical) = x.val

struct Superoxide <: Anion val end
struct Succinate <: Anion val end
struct Sodium <: Cation1 val end
struct Potassium <: Cation1 val end
struct Proton <: Cation1 val end
struct Calcium <: Cation2 val end
struct Magnesium <: Cation2 val end
struct ATP <: Adenylate val end
struct ADP <: Adenylate val end
struct AMP <: Adenylate val end
struct HATP <: Adenylate val end
struct HADP<: Adenylate  val end
struct MgATP <: Adenylate val end
struct MgADP <: Adenylate val end
struct PO4 <: Phosphate val end
struct HPO4 <: Phosphate val end
struct H2PO4 <: Phosphate val end

"""Ion valence"""
valence(x::Chemical) = 0
valence(x::Cation1) = +1
valence(x::Cation2) = +2
valence(x::Superoxide) = -1
valence(x::ATP) = -4
valence(x::ADP) = -3
valence(x::AMP) = -2
valence(x::PO4) = -3
valence(x::HATP) = valence(ATP(x)) + valence(Proton(x))
valence(x::HADP) = valence(ADP(x)) + valence(Proton(x))
valence(x::MgATP) = valence(ATP(x)) + valence(Magnesium(x))
valence(x::MgADP) = valence(ADP(x)) + valence(Magnesium(x))
valence(x::HPO4) = valence(PO4(x)) + valence(Proton(x))
valence(x::H2PO4) = valence(PO4(x)) + 2 * valence(Proton(x))

####################
# Binding constants
###################
#
# Abstract typefor binding polynomials
abstract type AffinityConstants end

"Dissociation constants from Guynn and Veech, 1973. T = 310K, Ionic strength = 0.2 M"
struct Guynn1973 <: AffinityConstants end

"Dissociation constants from Malyala et al. 2019. T = 310K, ionic strength = 0.17M"
struct Malyala2019 end

"Dissociation constants from Wei et al. 2011"
struct Wei2011 end

"Dissociation constants from O'Hara et al. 2011"
struct OHara2011 end

"""
Affinity constant for H + A = HA
"""
affinity(x::Chemical, y::Chemical, tag::AffinityConstants = Guynn1973()) = 0.0  # defaults to 0
affinity(y::Anion, x::Cation, tag::AffinityConstants) = affinity(x, y, tag)     # Order of cation and anion could be reversed
affinity(x::Proton, y::HPO4, tag::Guynn1973) = inv(1.76E-7M)
affinity(x::Proton, y::ATP, tag::Guynn1973) = inv(1.08E-7M)
affinity(x::Proton, y::ADP, tag::Guynn1973) = inv(1.20E-7M)
affinity(x::Proton, y::AMP, tag::Guynn1973) = inv(3.24E-7M)
affinity(x::Magnesium, y::ATP, tag::Guynn1973) = 13900/M
affinity(x::Magnesium, y::HATP, tag::Guynn1973) = 35.5/M
affinity(x::Magnesium, y::ADP, tag::Guynn1973) = 1320/M
affinity(x::Magnesium, y::HADP, tag::Guynn1973) = 32.4/M
affinity(x::Magnesium, y::AMP, tag::Guynn1973) = 60.1/M
affinity(x::Magnesium, y::HPO4, tag::Guynn1973) = 90.31/M
affinity(x::Proton, y::ATP, tag::Malyala2019) = inv(2.95E-7M)
affinity(x::Proton, y::ADP, tag::Malyala2019) = inv(4.37E-7M)
affinity(x::Proton, y::HPO4, tag::Malyala2019) = inv(2.19E-7M)
affinity(x::Sodium, y::ATP, tag::Malyala2019) = inv(72.4mM)
affinity(x::Sodium, y::ADP, tag::Malyala2019) = inv(105mM)
affinity(x::Sodium, y::HPO4, tag::Malyala2019) = inv(245mM)
affinity(x::Potassium, y::ATP, tag::Malyala2019) = inv(102mM)
affinity(x::Potassium, y::ADP, tag::Malyala2019) = inv(135mM)
affinity(x::Potassium, y::HPO4, tag::Malyala2019) = inv(389mM)
affinity(x::Magnesium, y::ATP, tag::Malyala2019) = inv(120μM)
affinity(x::Magnesium, y::ADP, tag::Malyala2019) = inv(955μM)
affinity(x::Magnesium, y::HPO4, tag::Malyala2019) = inv(33.9mM)
affinity(x::Calcium, y::ATP, tag::Malyala2019) = inv(129μM)
affinity(x::Calcium, y::ADP, tag::Malyala2019) = inv(1.51mM)
affinity(x::Calcium, y::HPO4, tag::Malyala2019) = inv(5.37mM)
affinity(x::Proton, y::ATP, tag::Wei2011) = inv(3.31E-7M)
affinity(x::Proton, y::ADP, tag::Wei2011) = inv(4.17E-7M)
affinity(x::Proton, y::HPO4, tag::Wei2011) = inv(1.70E-7M)
affinity(x::Proton, y::Succinate, tag::Wei2011) = inv(6.31E-6M)
affinity(x::Magnesium, y::ATP, tag::Wei2011) = inv(64.6μM)
affinity(x::Magnesium, y::ADP, tag::Wei2011) = inv(562μM)
affinity(x::Proton, y::HPO4, tag::OHara2011) = inv(1.76E-7M)
affinity(x::Sodium, y::HPO4, tag::OHara2011) = inv(224mM)
affinity(x::Potassium, y::HPO4, tag::OHara2011) = inv(292mM)

"Binding weight"
binding(x::Cation, y::Anion, tag::AffinityConstants = Guynn1973()) = qty(x) * affinity(x, y, tag)
binding(y::Anion, x::Cation, tag::AffinityConstants = Guynn1973()) = binding(x, y, tag)

"Binding polynomial, sum of binding weights"
poly(x::Cation, y::Anion, tag::AffinityConstants = Guynn1973()) = 1 + binding(x, y, tag)
poly(x1::Cation, x2::Cation, y::Anion, tag::AffinityConstants = Guynn1973()) = ploy(x1, y, tag) + binding(x2, y, tag)
poly(x1::Cation, x2::Cation, x3::Cation, y::Anion, tag::AffinityConstants = Guynn1973()) = ploy(x1, y, tag) + binding(x2, y, tag) + binding(x3, y, tag)

"Adenine nucleotide components"
function factions(h::Proton, mg::Magnesium, a::Adenylate, tag::AffinityConstants = Guynn1973())
    fH = binding(h, a, tag)
    fMg = binding(mg, a, tag)
    p = 1 + fH + fMg
    axpn = qty(a) / p
    return (axpn = axpn, haxp = fH * axpn, mgaxp = fMg * axpn, p = p)
end

"Phosphate components"
function factions(h::Proton, pi::Phosphate,tag::AffinityConstants = Guynn1973())
    p = poly(h, pi, tag)
    hpo4 = qty(pi) / p
    h2po4 = qty(pi) - hpo4
    return (hpo4 = hpo4, h2po4 = h2po4, p = p)
end
