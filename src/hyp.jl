# exported hyperbolic functions

over_sch(::Type{Float64}) = 710.0
over_sch(::Type{Float32}) = 89f0

"""
    sinh(x)

Compute hyperbolic sine of `x`.
"""
function sinh(x::T) where {T <: IEEEFloat}
    u = abs(x)
    d = expk2(Double(u))
    d = dsub(d, ddrec(d))
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
function cosh(x::T) where {T <: IEEEFloat}
    u = abs(x)
    d = expk2(Double(u))
    d = dadd(d, ddrec(d))
    u = T(d) * T(0.5)
    u = abs(x) > over_sch(T) ? T(Inf) : u
    u = isnan(u) ? T(Inf) : u
    u = isnan(x) ? T(NaN) : u
    return u
end



over_th(::Type{Float64}) = 18.714973875
over_th(::Type{Float32}) = 8.664339742f0

"""
    tanh(x)

Compute hyperbolic tangent of `x`.
"""
function tanh(x::T) where {T <: IEEEFloat}
    u = abs(x)
    d = expk2(Double(u))
    e = ddrec(d)
    d = ddiv(dsub(d, e), dadd(d, e))
    u = T(d)
    u = abs(x) > over_th(T) ? T(1) : u
    u = isnan(u) ? T(1) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u
    return u
end



"""
    asinh(x)

Compute the inverse hyperbolic sine of `x`.
"""
function asinh(x::T) where {T <: IEEEFloat}
    u = abs(x)
    d = logk2(dadd(dsqrt(dadd2(dsqu(u), T(1))), u))
    u = T(d)
    u = isinf(x) || isnan(u) ? T(Inf) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u
    return u
end



"""
    acosh(x)

Compute the inverse hyperbolic cosine of `x`.
"""
function acosh(x::T) where {T <: IEEEFloat}
    d = logk2(dadd2(dsqrt(dsub2(dsqu(x), T(1))), x))
    u = T(d)
    u = isinf(x) || isnan(u) ? T(Inf) : u
    u = x == T(1) ? T(0) : u
    u = x < T(1) ? T(NaN) : u
    u = isnan(x) ? T(NaN) : u
    return u
end



"""
    atanh(x)

Compute the inverse hyperbolic tangent of `x`.
"""
function atanh(x::T) where {T <: IEEEFloat}
    u = abs(x)
    d = logk2(ddiv(dadd2(T(1), u), dsub2(T(1), u)))
    u = u > T(1) ? T(NaN) : (u == T(1) ? T(Inf) : T(d) * T(0.5))
    u = isinf(x) || isnan(u) ? T(NaN) : u
    u = flipsign(u, x)
    u = isnan(x) ? T(NaN) : u
    return u
end
