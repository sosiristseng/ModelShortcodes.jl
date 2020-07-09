using ModelShortcodes
using Test

using Base.Iterators: product
import ModelShortcodes: Calcium,
                        Sodium,
                        Potassium,
                        Proton,
                        Superoxide,
                        Magnesium,
                        ATP,
                        ADP,
                        AMP,
                        HATP,
                        HADP,
                        MgATP,
                        MgADP,
                        Hydroxide,
                        Water,
                        Phosphate,
                        HPO4,
                        H2PO4

const xRange = 0.0:0.1:10.0
const kRange = 0.1:0.1:10.0
const nRange = 0.1:0.1:10.0

@testset "Valence" begin
    @test valence(Sodium) == +1
    @test valence(Potassium) == +1
    @test valence(Proton) == +1
    @test valence(Calcium) == +2
    @test valence(Superoxide) == -1
    @test valence(Hydroxide) == -1
    @test valence(Magnesium) == +2
    @test valence(ATP) == -4
    @test valence(ADP) == -3
    @test valence(AMP) == -2
    @test valence(Phosphate) == -3
    @test valence(HATP) == -3
    @test valence(HADP) == -2
    @test valence(MgATP) == -2
    @test valence(MgADP) == -1
    @test valence(HPO4) == -2
    @test valence(H2PO4) == -1
    @test valence(Water) == 0
    @test valence(Chemical) == 0
end

@testset "Functional" begin
    @test typeof(p_one(1.0)) === typeof(1.0)
    @test typeof(p_one(1)) === typeof(1)
    @test typeof(one_m(1.0)) === typeof(1.0)
    @test typeof(one_m(1)) === typeof(1)

    @test all(product(xRange, kRange)) do (x ,k) mm(x, k) ≈ x / (x + k) end
    @test all(product(xRange, kRange)) do (x ,k) mmr(x, k) ≈ mm(k, x) end
    @test all(product(xRange, kRange, nRange)) do (x, k, n) hill(x, k, n) ≈ mm(x^n, k^n) end
    @test all(product(xRange, kRange, nRange)) do (x, k, n) hillr(x, k, n) ≈ mm(k^n, x^n) end
    @test all(map(x -> expit(x) ≈ 1 / (1 + exp(-x)), -5.0:0.1:5.0))
    @test all(map(x -> exprel(x) ≈ x / (exp(x) - 1), kRange))
    @test all(map(x -> exprel(x) ≈ x / (exp(x) - 1), -kRange))
    @test exprel(0.0) ≈ 1.0
    @test all(map(x -> sqrt_s(x) ≈ sqrt(x), xRange))
    @test all(map(x -> sqrt_s(x) ≈ -sqrt(-x), -xRange))
    @test all(product(xRange, nRange)) do (x, n) pow_s(x, n) ≈ x^n end
    @test all(product(-xRange, nRange)) do (x, n) pow_s(x, n) ≈ -((-x)^n) end
    @test all(product(xRange, kRange, nRange)) do (x, k, n) hill_s(x, k, n) ≈ hill(x, k, n) end
    @test all(product(-xRange, kRange, nRange)) do (x, k, n) hill_s(x, k, n) ≈ -hill(-x, k, n) end

    iStim = i_stim.(0:20; period=10, duty=1, strength=2)
    @test (iStim[1], iStim[2], iStim[4], iStim[10], iStim[11]) == (2, 0, 0, 0, 2)
end


@testset "Nernst equation" begin
    @test nernst(1.0, 1.0, Calcium) ≈ zero(1.0)

    # https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html
    @test nernst(2.0, 70E-6, Calcium) ≈ 0.137 rtol=0.01  # Calcium
    @test nernst(4.0, 150.0) ≈ -0.09681 rtol=0.01        # Potassium
    @test pmf(0.0, 1E-4, 1E-4) ≈ 0.0                     # Complete uncoupled
    @test pmf(0.15, 1E-4, 5*1E-5) ≈ 0.169 rtol=0.01      # Physiological PMF
    @test pmf(0.15, -0.3) ≈ 0.168 rtol=0.01              # Physiological PMF
    @test ghk(1.0, Calcium, 0.0, 1.0, 1.0) ≈ 0.0
    @test ghk(1E-6, Calcium, -0.050, 70E-6, 2.0) ≈ -1.48 rtol=0.01
end
