@testset "Functional" begin

    # Make sure the Type does not change
    @inferred p_one(1.0)
    @inferred p_one(1)
    @inferred one_m(1.0)
    @inferred one_m(1)

    @testset "Hill functions" for x in 0.0:0.5:5.0,  k in 0.1:0.1:5.0
        @test mm(x, k) ≈ hill(x, k)
        @test mmr(x, k) ≈ hillr(x, k)
        for n in 0.1:0.5:5.0
            a = hill(x, k, n)
            b = hillr(x, k, n)
            @test a + b ≈ 1.0
            @test hill(x+0.1, k, n) >= a
            @test hill(x, k+0.1, n) <= a
        end
    end

    @test expit(-1e9) ≈ 0.0
    @test expit(0.0) ≈ 0.5
    @test expit(1e9) ≈ 1.0
    @test exprel(-1e9) ≈ 1e9 atol = 1e-5
    @test exprel(0.0) ≈ 1.0
    @test exprel(1e9) ≈ 0.0 atol = 1e-5

    @testset "Signed functions" for x in -0.1:-0.5:-6.0
        @test sqrt_s(x) ≈ - sqrt(-x)
        for n in 0.1:1.0:10.0
            @test pow_s(x, n) ≈ - ((-x)^n)
            for k in -10.0:1.0:10.0
                @test hill_s(x, k, n) ≈ hillr_s(k, x, n)
                @test hill_s(x, k, n) + hillr_s(x, k, n) ≈ sign(x) * sign(k)
            end
        end
    end

    @test stimulate.(0:20; period=10, duty=1, strength=2)[[1, 2, 4, 10, 11]] == [2, 0, 0, 0, 2]
end

@testset "Nernst and GHK" begin
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
