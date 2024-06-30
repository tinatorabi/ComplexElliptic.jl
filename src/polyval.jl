"""
    polyval(p::Vector{T}, x::T) where T <: Number

Evaluate a polynomial at a given point `x`. The coefficients of the polynomial are provided in the vector `p`.

# Arguments
- `p::Vector{T}`: A vector of coefficients.
- `x::T`: The point at which the polynomial is evaluated.

# Returns
- The value of the polynomial at `x`.
"""
function polyval(p::Vector{T}, x::T) where T <: Number
    result = zero(T)
    for i in 1:length(p)
        result += p[i] * x^(length(p) - i)
    end
    return result
end
