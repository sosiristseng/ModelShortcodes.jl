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
                        H2PO4,
                        HydrogenPeroxide

@testset "Chemical" begin
    @test valence(Sodium) == 1
    @test valence(Potassium) == 1
    @test valence(Proton) == 1
    @test valence(Calcium) == 2
    @test valence(Magnesium) == 2
    @test valence(Superoxide) == -1
    @test valence(Hydroxide) == -1
    @test valence(ATP) == -4
    @test valence(ADP) == -3
    @test valence(AMP) == -2
    @test valence(HATP) == -3
    @test valence(HADP) == -2
    @test valence(MgATP) == -2
    @test valence(MgADP) == -1
    @test valence(HPO4) == -2
    @test valence(H2PO4) == -1
    @test valence(Water) == 0
    @test valence(HydrogenPeroxide) == 0
end

@testset "Functional" begin

    for x in  (1.0, 1)
        @test typeof(p_one(x)) === typeof(x)
        @test typeof(one_m(x)) === typeof(x)
    end

    for x in (0, 1, 1e9)
        @test mm(x, 1) + mmr(x, 1) ≈ 1.0
        @test hill(x, 1, 4) + hillr(x, 1, 4) ≈ 1.0
    end

    @test expit(-1e9) ≈ 0.0
    @test expit(0.0) ≈ 0.5
    @test expit(1e9) ≈ 1.0
    @test exprel(-1e9) ≈ 1e9 atol = 1e-5
    @test exprel(0.0) ≈ 1.0
    @test exprel(1e9) ≈ 0.0 atol = 1e-5

    for x in (-1.0, 0.0, 1.0)
        @test sqrt_s(x) ≈ sign(x) * sqrt(abs(x))
        @test pow_s(x, 2) ≈ sign(x) * pow_s(abs(x), 2)
        @test hill_s(x, 1.0, 2.0) ≈ sign(x) * hill(abs(x), 1.0, 2.0)
    end

    @test i_stim.(0:20; period=10, duty=1, strength=2)[[1, 2, 4, 10, 11]] == [2, 0, 0, 0, 2]
end

@testset "Nernst equation" begin
    # Zero conditions
    @test nernst(1.0, 1.0, Calcium) ≈ zero(1.0)
    @test pmf(0.0, 1E-4, 1E-4) ≈ 0.0  # Complete uncoupled
    @test ghk(1.0, Calcium, 0.0, 1.0, 1.0) ≈ 0.0

    # Non-zero conditions
    # https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html
    @test nernst(2.0, 70E-6, Calcium) ≈ 0.137 rtol=0.01  # Calcium
    @test nernst(4.0, 150.0) ≈ -0.09681 rtol=0.01        # Potassium
    @test pmf(0.15, 1E-4, 5*1E-5) ≈ 0.169 rtol=0.01      # Physiological PMF
    @test pmf(0.15, -0.3) ≈ 0.168 rtol=0.01              # Physiological PMF
    @test ghk(1E-6, Calcium, -0.050, 70E-6, 2.0) ≈ -1.48 rtol=0.01
end


@testset "Binding" begin
    @test power_to_conc(7) ≈ 1E-4
    @test conc_to_power(1E-4) ≈ 7

    @test binding(Water, Water) |> isinf

    cations = (Proton, Magnesium, Sodium, Potassium, Calcium)
    anions = (Hydroxide, HPO4, ATP, ADP, AMP)

    # should be no binding bewteen ions and water
    @test binding(Water, Water) |> isinf
    for x in (cations..., anions...)
        @test binding(Water, x) |> isinf
        @test binding(x, Water) |> isinf
    end

    for cation in cations, anion in anions
        # The same chemical
        @test binding(cation, cation) |> isinf
        @test binding(anion, anion) |> isinf
        # Cations and anions should be switchable
        @test binding(cation, anion) ≈ binding(anion, cation)
    end

    h = power_to_conc(7)
    @test binding_fraction(h, HPO4, Proton) != 0.0
    @test binding_fraction(h, HPO4, Water) ≈ 0.0
    @test poly(h, HPO4, Water) ≈ 1.0
    @test poly(h, HPO4, Proton) > 1.0
end
