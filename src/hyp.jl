# exported hyperbolic functions

over_sch(::Type{Float64}) = 710.0
over_sch(::Type{Float32}) = 89f0

"""
    sinh(x)

Compute hyperbolic sine of `x`.
"""
function sinh(x::T) where {T<:Union{Float32,Float64}}
    u = abs(x)
    d = expk2(Double(u))
    d = dsub(d, drec(d))
    u = T(d) * T(0.5)
    u = abs(x) > over_sch(T) ? T(Inf) : u
    u = isnan(u) ? T(Inf) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u
    return u
end



"""
    cosh(x)

Compute hyperbolic cosine of `x`.
"""
function cosh(x::T) where {T<:Union{Float32,Float64}}
    u = abs(x)
    d = expk2(Double(u))
    d = dadd(d, drec(d))
    u = T(d) * T(0.5)
    u = abs(x) > over_sch(T) ? T(Inf) : u
    u = isnan(u) ? T(Inf) : u
    u = isnan(x) ? T(NaN) : u
    return u
end



over_th(::Type{Float64}) = 18.714973875
over_th(::Type{Float32}) = 18.714973875f0

"""
    tanh(x)

Compute hyperbolic tangent of `x`.
"""
function tanh(x::T) where {T<:Union{Float32,Float64}}
    u = abs(x)
    d = expk2(Double(u))
    e = drec(d)
    d = ddiv(dsub(d, e), dadd(d, e))
    u = T(d)
    u = abs(x) > over_th(T) ? T(1.0) : u
    u = isnan(u) ? T(1) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u
    return u
end



"""
    asinh(x)

Compute the inverse hyperbolic sine of `x`.
"""
function asinh(x::T) where {T<:Union{Float32,Float64}}
    y = abs(x)

    d = y > 1 ? drec(x) : Double(y, T(0.0))
    d = dsqrt(dadd2(dsqu(d), T(1.0)))
    d = y > 1 ? dmul(d, y) : d

    d = logk2(dnormalize(dadd(d, x)))
    y = T(d)

    y = (abs(x) > SQRT_MAX(T) || isnan(y)) ? flipsign(T(Inf), x) : y
    y = isnan(x) ? T(NaN) : y
    y = isnegzero(x) ? T(-0.0) : y
    
    return y
end



"""
    acosh(x)

Compute the inverse hyperbolic cosine of `x`.
"""
function acosh(x::T) where {T<:Union{Float32,Float64}}
    d = logk2(dadd2(dmul(dsqrt(dadd2(x, T(1.0))), dsqrt(dsub2(x, T(1.0)))), x))
    y = T(d)

    y = (x > SQRT_MAX(T) || isnan(y)) ? T(Inf) : y
    y = x == T(1.0) ? T(0.0) : y
    y = x < T(1.0) ? T(NaN) : y
    y = isnan(x) ? T(NaN) : y

    return y
end



"""
    atanh(x)

Compute the inverse hyperbolic tangent of `x`.
"""
function atanh(x::T) where {T<:Union{Float32,Float64}}
    u = abs(x)
    d = logk2(ddiv(dadd2(T(1.0), u), dsub2(T(1.0), u)))
    u = u > T(1.0) ? T(NaN) : (u == T(1.0) ? T(Inf) : T(d) * T(0.5))

    u = isinf(x) || isnan(u) ? T(NaN) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u

    return u
end
