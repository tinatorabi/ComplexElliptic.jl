
module ComplexElliptic

"""
MATLAB Code:
Author: Toby Driscoll
Year: 1999


Julia Rewrite:
Author: Tina Torabi
Year: 2023
"""


using SpecialFunctions, Elliptic

function polyval(p::Vector{T}, x::T) where T <: Real
    result = zero(T) 
    for i in 1:length(p)
        result += p[i] * x^(length(p) - i)
    end
    return result
end

function ellipjc(u, L, flag=false)
    u = complex(u)
    if !flag
	L_= exp(-2 * π * L)
        K, Kp = ellipke(L_)[1], ellipke(1-L_)[1]
        high = imag(u) > Kp / 2
        if high
            u = im * Kp - u
        end
        m = exp(-2 * pi * L)
    else
        high = false
        m = L
    end

    if m < 6eps(Float64)
        sinu = sin(u)
        cosu = cos(u)
        sn = sinu + m / 4 * (sinu .* cosu - u) .* cosu
        cn = cosu - m / 4 * (sinu .* cosu - u) .* sinu
        dn = 1 - m / 4 * (sinu .^ 2 - cosu .^ 2)
    else
        kappa = m > 1e-3 ? (1 - sqrt(1 - m)) / (1 + sqrt(1 - m)) : polyval([132.0, 42.0, 14.0, 5.0, 2.0, 1.0, 0.0], m / 4.0)
        mu = kappa^2
        v = u / (1 + kappa)
        sn1, cn1, dn1 = ellipjc(v, mu, true)
        
        denom = 1 + kappa .* sn1 .^ 2
        sn = (1 + kappa) .* sn1 ./ denom
        cn = cn1 .* dn1 ./ denom
        dn = (1 - kappa .* sn1 .^ 2) ./ denom
    end

    if high
        sn = -1 / (sqrt(m) * sn)
        cn = im * dn / (sqrt(m) * sn)
        dn = im * cn / sn
    end

    return sn, cn, dn
end


end # module ComplexElliptic