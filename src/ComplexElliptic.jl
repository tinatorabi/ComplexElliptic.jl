module ComplexElliptic

"""

Author: Tina Torabi
Year: 2023
"""

function polyval(p::Vector{T}, x::T) where T <: Number
    result = zero(T) 
    for i in 1:length(p)
        result += p[i] * x^(length(p) - i)
    end
    return result
end



 function ellipkkp(L)
    if abs(L) > 10
        K = π / 2
        Kp = π * L + log(4)
        return K, Kp
    end

    m = exp(-2 * π * L)  
    a0 = 1.0
    b0 = sqrt(1 - m)
    a1=0
    s0 = m
    i1 = 0
    mm = 1.0
    while abs(mm) > eps()
        a1 = (a0 + b0) / 2
        b1 = sqrt(a0 * b0)
        c1 = (a0 - b0) / 2
        i1 += 1
        w1 = 2^i1 * c1^2
        mm = w1[argmax(abs.(w1))]  
        s0 += w1
        a0 = a1
        b0 = b1
    end
    K_ = π / (2 * a1)

    a0 = 1.0
    b0 = sqrt(m)
    a1 = 0
    s0 = 1 - m
    i1 = 0
    mm = 1.0
    while abs(mm) > eps()
        a1 = (a0 + b0) / 2
        b1 = sqrt(a0 * b0)
        c1 = (a0 - b0) / 2
        i1 += 1
        w1 = 2^i1 * c1^2
        mm = w1[argmax(abs.(w1))]  
        s0 += w1
        a0 = a1
        b0 = b1
    end
    Kp = π / (2 * a1)

    if m == 0
        Kp = Inf
    end

    return K_, Kp
end



function ellipjc(u, L; flag=false)
    u = complex(u)
    if !flag
        KK, Kp = ellipkkp(L)
        high = imag.(u) .> real(Kp) / 2
        u[high] = im * Kp .- u[high]
        m = exp(-2 * pi * L)
    else
        high = falses(size(u))  
        m = L
    end

    if abs(m) < 6eps(Float64)
        sinu = sin.(u)
        cosu = cos.(u)
        sn = sinu + m / 4 * (sinu .* cosu - u) .* cosu
        cn = cosu - m / 4 * (sinu .* cosu - u) .* sinu
        dn = 1 .+ m / 4 .* (cosu .^ 2 - sinu .^ 2 .- 1)
    else
        if abs(m) > 1e-3
            kappa = (1 - sqrt(Complex(1 - m))) / (1 + sqrt(Complex(1 - m)))
        else
            kappa = polyval(Complex{Float64}.([132.0, 42.0, 14.0, 5.0, 2.0, 1.0, 0.0]), Complex{Float64}(m / 4.0))
        end
        mu = kappa ^ 2
        v = u ./ (1 + kappa)
        sn1, cn1, dn1 = ellipjc(v, mu, flag=true)

        denom = 1 .+ kappa .* sn1 .^ 2
        sn = (1 .+ kappa) .* sn1 ./ denom
        cn = cn1 .* dn1 ./ denom
        dn = (1 .- kappa .* sn1 .^ 2) ./ denom
    end

    if any(high)
        snh = sn[high]
        cnh = cn[high]
        dnh = dn[high]
        sn[high] = -1 ./ (sqrt(m) * snh)
        cn[high] = im * dnh ./ (sqrt(m) * snh)
        dn[high] = im * cnh ./ snh
    end

    return sn, cn, dn
end

end # module ComplexElliptic
