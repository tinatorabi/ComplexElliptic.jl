"""
    ellipjc(u::Complex, L::Real)

Compute the Jacobi elliptic functions `sn`, `cn`, and `dn` for a complex argument `u` and
parameter `m = exp(-2*pi*L)`, where `0 < L < Inf`. The relationship `m = k^2` applies,
where `k` is the elliptic modulus.

# Arguments
- `u::Complex`: A complex number or an array of complex numbers. Each element should
  satisfy the condition `|Re(u)| < K` for real part and `0 < Im(u) < Kp` for imaginary
  part, where `[K, Kp] = ellipkkp(L)`.
- `L::Real`: A real scalar defining the modulus squared parameter `m`.

# Returns
- `(sn, cn, dn)`: Corresponding to the Jacobi elliptic functions evaluated at each element of `u`.

# Examples
```julia
L = 0.5
u = 1.0 + 0.5im
sn, cn, dn = ellipjc(u, L)
"""
function ellipjc(u, L; flag=false)
    if !flag
        _, Kp = ellipkkp(L)
        high = [(imag(x) > real(Kp) / 2) for x in u]
        if any(high)
            u = [(imag(x) > real(Kp) / 2) ? im * Kp .- x : x for x in u]
        end
        m = exp(-2 * pi * L)
    else
        high = falses(size(u))  
        m = L
    end
    if abs(m) < 6eps(Float64)
        sinu = [sin(x) for x in u]
        cosu = [cos(x) for x in u]
        sn = sinu + m / 4 * (sinu .* cosu - u) .* cosu
        cn = cosu - m / 4 * (sinu .* cosu - u) .* sinu
        dn = 1 .+ m / 4 .* (cosu .^ 2 - sinu .^ 2 .- 1)
    else
        if abs(m) > 1e-3
            kappa = (1 - sqrt(Complex(1 - m))) / (1 + sqrt(Complex(1 - m)))
        else
            kappa = ComplexElliptic.polyval(Complex{Float64}.([132.0, 42.0, 14.0, 5.0, 2.0, 1.0, 0.0]), Complex{Float64}(m / 4.0))
        end
        sn1, cn1, dn1 = ellipjc(u / (1 + kappa), kappa ^ 2, flag=true)

        denom = 1 .+ kappa .* sn1 .^ 2
        sn = (1 .+ kappa) .* sn1 ./ denom
        cn = cn1 .* dn1 ./ denom
        dn = (1 .- kappa .* sn1 .^ 2) ./ denom
    end

    if any(high)
        snh = Zygote.Buffer(sn)
        cnh = Zygote.Buffer(cn)
        dnh = Zygote.Buffer(dn)
        snh[:] = sn[:]
        cnh[:] = cn[:]
        dnh[:] = dn[:]
        snh[high] = -1 ./ (sqrt(m) * sn[high])
        cnh[high] = im * dn[high] ./ (sqrt(m) * sn[high])
        dnh[high] = im * cn[high] ./ sn[high]
    end

    return copy(sn), copy(cn), copy(dn)
end
