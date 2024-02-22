using Test
using ComplexElliptic  

# Test polyval function
@testset "polyval Function" begin
    p = [1.0, 0.0, -1.0] # Represents the polynomial x^2 - 1
    x = 0.0
    result = ComplexElliptic.polyval(p, x)
    @test result ≈ -1.0
end


# Test ellipjc function 
# comparing to values obtained from @Toby Driscoll

@testset "ellipjc Function" begin
    u = 0.5
    L = 1.0
    sn, cn, dn = ComplexElliptic.ellipjc(u, L)
    # Verify the types of the returned values
    @test typeof(sn) == Complex{Float64}
    @test typeof(cn) == Complex{Float64}
    @test typeof(dn) == Complex{Float64}
    @test round(abs(sn[1]), digits=4) ≈ 0.4794
    @test round(abs(cn[1]), digits=4) ≈ 0.8776
    @test round(abs(dn[1]), digits=4) ≈ 0.9998
end
