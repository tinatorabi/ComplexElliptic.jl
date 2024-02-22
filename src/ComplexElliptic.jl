module ComplexElliptic

"""

Author: Tina Torabi
Year: 2023
"""

function ellipkkp(L)
    # Complete elliptic integral of the first kind, with complement.
    if L > 10
        K = π / 2
        Kp = π * L + log(4)
        return K, Kp
    end

    m = exp(-2 * π * L)
    a0 = 1.0
    b0 = sqrt(1 - m)
    s0 = m
    i1 = 0
    mm = 1.0
    while mm > eps()
        a1 = (a0 + b0) / 2
        b1 = sqrt(a0 * b0)
        c1 = (a0 - b0) / 2
        i1 += 1
        w1 = 2^i1 * c1^2
        mm = maximum(w1)
        s0 += w1
        a0 = a1
        b0 = b1
    end
    K = π / (2 * a1)

    if m == 1
        K = Inf
    end

    if nargout() > 1
        a0 = 1.0
        b0 = sqrt(m)
        a1=0
        s0 = 1 - m
        i1 = 0
        mm = 1.0
        while mm > eps()
            a1 = (a0 + b0) / 2
            b1 = sqrt(a0 * b0)
            c1 = (a0 - b0) / 2
            i1 += 1
            w1 = 2^i1 * c1^2
            mm = maximum(w1)
            s0 += w1
            a0 = a1
            b0 = b1
        end
        Kp = π / (2 * a1)
        if m == 0
            Kp = Inf
        end
    else
        Kp = nothing
    end

    return K, Kp
end


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
        K, Kp = ellipkkp(L_)[1], ellipkkp(1-L_)[1]
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
