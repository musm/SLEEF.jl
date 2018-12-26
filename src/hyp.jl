# exported hyperbolic functions

over_sch(::Type{Float64}) = 710.0
over_sch(::Type{Float32}) = 89f0

"""
    sinh(x)

Compute hyperbolic sine of `x`.
"""
function sinh(x::V) where {V <: FloatType}
    T = eltype(x)
    u = abs(x)
    d = expk2(Double(u))
    d = dsub(d, drec(d))
    u = V(d) * T(0.5)
    u = vifelse(abs(x) > over_sch(T), T(Inf), u)
    u = vifelse(isnan(u), T(Inf), u)
    u = flipsign(u, x)
    u = vifelse(isnan(x), T(NaN), u)
    return u
end



"""
    cosh(x)

Compute hyperbolic cosine of `x`.
"""
function cosh(x::V) where {V <: FloatType}
    T = eltype(x)
    u = abs(x)
    d = expk2(Double(u))
    d = dadd(d, drec(d))
    u = V(d) * T(0.5)
    u = vifelse(abs(x) > over_sch(T), T(Inf), u)
    u = vifelse(isnan(u), T(Inf), u)
    u = vifelse(isnan(x), T(NaN), u)
    return u
end



over_th(::Type{Float64}) = 18.714973875
over_th(::Type{Float32}) = 18.714973875f0

"""
    tanh(x)

Compute hyperbolic tangent of `x`.
"""
function tanh(x::V) where {V <: FloatType}
    T = eltype(x)
    u = abs(x)
    d = expk2(Double(u))
    e = drec(d)
    d = ddiv(dsub(d, e), dadd(d, e))
    u = V(d)
    u = vifelse(abs(x) > over_th(T), T(1.0), u)
    u = vifelse(isnan(u), T(1), u)
    u = flipsign(u, x)
    u = vifelse(isnan(x), T(NaN), u)
    return u
end



"""
    asinh(x)

Compute the inverse hyperbolic sine of `x`.
"""
function asinh(x::V) where {V <: FloatType}
    T = eltype(x)
    y = abs(x)

    yg1 = y > 1
    d = vifelse(yg1, drec(x), Double(y, V(0.0)))
    d = dsqrt(dadd2(dsqu(d), T(1.0)))
    d = vifelse(yg1, dmul(d, y), d)

    d = logk2(dnormalize(dadd(d, x)))
    y = V(d)

    y = vifelse(((abs(x) > SQRT_MAX(T)) | isnan(y)), flipsign(T(Inf), x), y)
    y = vifelse(isnan(x), T(NaN), y)
    y = vifelse(isnegzero(x), T(-0.0), y)

    return y
end



"""
    acosh(x)

Compute the inverse hyperbolic cosine of `x`.
"""
function acosh(x::V) where {V <: FloatType}
    T = eltype(x)
    d = logk2(dadd2(dmul(dsqrt(dadd2(x, T(1.0))), dsqrt(dsub2(x, T(1.0)))), x))
    y = V(d)

    y = vifelse(((x > SQRT_MAX(T)) | isnan(y)), T(Inf), y)
    y = vifelse(x == T(1.0), T(0.0), y)
    y = vifelse(x < T(1.0), T(NaN), y)
    y = vifelse(isnan(x), T(NaN), y)

    return y
end



"""
    atanh(x)

Compute the inverse hyperbolic tangent of `x`.
"""
function atanh(x::V) where {V <: FloatType}
    T = eltype(x)
    u = abs(x)
    d = logk2(ddiv(dadd2(T(1.0), u), dsub2(T(1.0), u)))
    u = vifelse(u > T(1.0), T(NaN), vifelse(u == T(1.0), T(Inf), V(d) * T(0.5)))

    u = vifelse(isinf(x) | isnan(u), T(NaN), u)
    u = flipsign(u, x)
    u = vifelse(isnan(x), T(NaN), u)

    return u
end
