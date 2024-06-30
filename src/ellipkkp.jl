"""
    ellipkkp(L) -> (K, Kp)

Calculate the complete elliptic integrals of the first kind (K) and its complement (Kp) using the parameter `L`.
The function evaluates these integrals at `M = exp(-2 * π * L)`, where `0 < L < ∞`.

# Arguments
- `L::Real`: A positive real number that parameterizes the modulus of the elliptic integral.

# Returns
- `(K, Kp)`: A tuple where `K` is the complete elliptic integral of the first kind for the modulus `M`, and `Kp` is the complete elliptic integral for the complementary modulus `1 - M`.

# Details
The function computes the complete elliptic integrals where:
- `M = exp(-2 * π * L)` serves as the elliptic modulus squared (`k^2`).

# Examples
```julia
K, Kp = ellipkkp(1.0)
"""
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
