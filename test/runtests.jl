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
                        Phosphate,
                        HPO4,
                        H2PO4

@testset "ModelShortcodes" begin

@testset "Chemical" begin
    x = zero(Float64)
    @test valence(Sodium(x)) == 1
    @test valence(Potassium(x)) == 1
    @test valence(Proton(x)) == 1
    @test valence(Calcium(x)) == 2
    @test valence(Magnesium(x)) == 2
    @test valence(Superoxide(x)) == -1
    @test valence(ATP(x)) == -4
    @test valence(ADP(x)) == -3
    @test valence(AMP(x)) == -2
    @test valence(HATP(x)) == -3
    @test valence(HADP(x)) == -2
    @test valence(MgATP(x)) == -2
    @test valence(MgADP(x)) == -1
    @test valence(HPO4(x)) == -2
    @test valence(H2PO4(x)) == -1
end

@testset "Functional" begin

    @inferred p_one(1.0)
    @inferred p_one(1)
    @inferred one_m(1.0)
    @inferred one_m(1)

    @testset "Hill functions" for x in 0.0:1.0:10.0
        @test mm(x) ≈ hill(x)
        @test mm(x) ≈ mm(x, 1.0)
        @test mmr(x) ≈ hillr(x)
        @test mmr(x) ≈ mmr(x, 1.0)
        for k in 0.1:1.0:10.0
            @test mm(x, k) ≈ hill(x, k, 1)
            @test mmr(x, k) ≈ hillr(x, k, 1)
            @test mm(x, k) ≈ mmr(k, x)
            @test mm(x, k) + mmr(x, k) ≈ 1.0
            for n in 0.1:1.0:10.0
                @test hill(x, k, n) ≈ hillr(k, x, n)
                @test hill(x, k, n) + hillr(x, k, n) ≈ 1.0
            end
        end
    end

    @test expit(-1e9) ≈ 0.0
    @test expit(0.0) ≈ 0.5
    @test expit(1e9) ≈ 1.0
    @test exprel(-1e9) ≈ 1e9 atol = 1e-5
    @test exprel(0.0) ≈ 1.0
    @test exprel(1e9) ≈ 0.0 atol = 1e-5

    @testset "Signed functions" for x in -0.1:-1.0:-10.0
        @test sqrt_s(x) ≈ - sqrt(-x)
        for n in 0.1:1.0:10.0
            @test pow_s(x, n) ≈ - ((-x)^n)
            for k in -10.0:1.0:10.0
                @test hill_s(x, k, n) ≈ hillr_s(k, x, n)
                @test hill_s(x, k, n) + hillr_s(x, k, n) ≈ sign(x) * sign(k)
            end
        end
    end

    @test i_stim.(0:20; period=10, duty=1, strength=2)[[1, 2, 4, 10, 11]] == [2, 0, 0, 0, 2]
end

@testset "Nernst equation" begin
    ca = Calcium
    # Zero conditions
    @test nernst(1.0, 1.0, 2) ≈ zero(1.0)
    @test pmf(0.0, 1E-4, 1E-4) ≈ 0.0  # Complete uncoupled
    @test ghk(1.0, 2, 0.0, 1.0, 1.0) ≈ 0.0

    # Non-zero conditions
    # https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html
    @test nernst(2.0, 70E-6, 2) ≈ 0.137 rtol=0.01  # Calcium
    @test nernst(4.0, 150.0) ≈ -0.09681 rtol=0.01        # Potassium
    @test pmf(0.15, 1E-4, 5*1E-5) ≈ 0.169 rtol=0.01      # Physiological PMF
    @test pmf(0.15, -0.3) ≈ 0.168 rtol=0.01              # Physiological PMF
    @test ghk(1E-6, 2, -0.050, 70E-6, 2.0) ≈ -1.48 rtol=0.01
end


@testset "Binding" begin
    @test power_to_conc(7) ≈ 1E-4
    @test conc_to_power(1E-4) ≈ 7
end

end # Test set ModelShortcodes