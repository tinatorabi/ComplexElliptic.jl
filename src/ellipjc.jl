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
    u = u isa Array ? u : [u]
    sn = similar(u, Complex{Float64})
    cn = similar(u, Complex{Float64})
    dn = similar(u, Complex{Float64})
    
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
        sn .= sinu + m / 4 * (sinu .* cosu - u) .* cosu
        cn .= cosu - m / 4 * (sinu .* cosu - u) .* sinu
        dn .= 1 .+ m / 4 .* (cosu .^ 2 - sinu .^ 2 .- 1)
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
        sn .= (1 .+ kappa) .* sn1 ./ denom
        cn .= cn1 .* dn1 ./ denom
        dn .= (1 .- kappa .* sn1 .^ 2) ./ denom
    end

    if any(high)
        snh = sn[high]
        cnh = cn[high]
        dnh = dn[high]
        sn[high] .= -1 ./ (sqrt(m) * snh)
        cn[high] .= im * dnh ./ (sqrt(m) * snh)
        dn[high] .= im * cnh ./ snh
    end

    return sn, cn, dn
end
