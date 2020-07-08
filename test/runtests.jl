using ModelShortcodes
using Test

@testset "ModelShortcodes.jl" begin
    @testset "Valences" begin
        import ModelShortcodes: Calcium, Sodium, Potassium, Proton, Superoxide, Magnesium, ATP, ADP, AMP, HATP, HADP, MgATP, MgADP, Hydroxide, Water, Phosphate, HPO4, H2PO4
        @test valence(Water) == 0
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
    end

    @testset "Functional" begin
        @test isa(p_one(1.0), typeof(1.0))
        @test isa(p_one(1), typeof(1))
        @test p_one(0.5) ≈ 1.5
        @test isa(one_m(1.0), typeof(1.0))
        @test isa(one_m(1), typeof(1))
        @test one_m(0.5) ≈ 0.5
        @test mm(0.0) ≈ 0.0
        @test mm(1.0) ≈ 0.5
        @test mm(3.0, 3.0) ≈ 0.5
        @test mm(10000.0, 1.0) > 0.99
        @test all(Base.Iterators.product(0.0:0.1:10.0, 0.1:0.1:10.0)) do (x ,k)
            mmr(x, k) + mm(x, k) ≈ 1.0
        end
        @test all(Base.Iterators.product(0.0:0.1:10.0, 0.1:0.1:10.0, 0.1:0.1:10.0)) do (x ,k, n)
            hill(x, k, n) + hillr(x, k, n) ≈ 1.0
        end
    end
end
