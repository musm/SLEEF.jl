
"""
  pow(x, y)

Exponentiation operator, returns `x` raised to the power `y`.
"""
function pow(x::T, y::T) where {T<:Union{Float32,Float64}}
    yi = unsafe_trunc(Int, y)
    yisint = yi == y
    yisodd = isodd(yi) && yisint
    
    result = expk(dmul(logk(abs(x)), y))

    result = isnan(result) ? T(Inf) : result
    result *= (x > 0 ? T(1.0) : (!yisint ? T(NaN) : (yisodd ? -T(1.0) : T(1.0))))

    efx = flipsign(abs(x) - 1, y)
    isinf(y) && (result = efx < 0 ? T(0.0) : (efx == 0 ? T(1.0) : T(Inf)))
    (isinf(x) || x == 0) && (result = (yisodd ? _sign(x) : T(1.0)) * ((x == 0 ? -y : y) < 0 ? T(0.0) : T(Inf)))
    (isnan(x) || isnan(y)) && (result = T(NaN))
    (y == 0 || x == 1) && (result = T(1.0))

    return result
end



let
global cbrt_fast
global cbrt

c6d = -0.640245898480692909870982
c5d = 2.96155103020039511818595
c4d = -5.73353060922947843636166
c3d = 6.03990368989458747961407
c2d = -3.85841935510444988821632
c1d = 2.2307275302496609725722

c6f = -0.601564466953277587890625f0
c5f =  2.8208892345428466796875f0
c4f = -5.532182216644287109375f0
c3f =  5.898262500762939453125f0
c2f = -3.8095417022705078125f0
c1f =  2.2241256237030029296875f0

global @inline cbrt_kernel(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d
global @inline cbrt_kernel(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f

"""
    cbrt_fast(x)

Return `x^{1/3}`.
"""
function cbrt_fast(d::T) where {T<:Union{Float32,Float64}}
    e  = ilogbk(abs(d)) + 1
    d  = ldexp2k(d, -e)
    r  = (e + 6144) % 3
    q  = r == 1 ? T(M2P13) : T(1)
    q  = r == 2 ? T(M2P23) : q
    q  = ldexp2k(q, (e + 6144) ÷ 3 - 2048)
    q  = flipsign(q, d)
    d  = abs(d)
    x  = cbrt_kernel(d)
    y  = x * x
    y  = y * y
    x -= (d * y - x) * T(1 / 3)
    y  = d * x * x
    y  = (y - T(2 / 3) * y * (y * x - 1)) * q
end


"""
    cbrt(x)

Return `x^{1/3}`. The prefix operator `∛` is equivalent to `cbrt`.
"""
function cbrt(d::T) where {T<:Union{Float32,Float64}}
    e  = ilogbk(abs(d)) + 1
    d  = ldexp2k(d, -e)
    r  = (e + 6144) % 3
    q2 = r == 1 ? MD2P13(T) : Double(T(1))
    q2 = r == 2 ? MD2P23(T) : q2
    q2 = flipsign(q2, d)
    d  = abs(d)
    x  = cbrt_kernel(d)
    y  = x * x
    y  = y * y
    x -= (d * y - x) * T(1 / 3)
    z  = x
    u  = dsqu(x)
    u  = dsqu(u)
    u  = dmul(u, d)
    u  = dsub(u, x)
    y  = T(u)
    y  = -T(2 / 3) * y * z
    v  = dadd(dsqu(z), y)
    v  = dmul(v, d)
    v  = dmul(v, q2)
    z  = ldexp2k(T(v), (e + 6144) ÷ 3 - 2048)
    isinf(d) && (z = flipsign(T(Inf), q2.hi))
    d == 0   && (z = flipsign(T(0), q2.hi))
    return z
end
end



"""
    hypot(x,y)

Compute the hypotenuse `\\sqrt{x^2+y^2}` avoiding overflow and underflow.
"""
function hypot(x::T, y::T) where {T<:IEEEFloat}
    x = abs(x)
    y = abs(y)
    if x < y
       x, y = y, x
    end
    r = (x == 0) ? y : y / x
    x * sqrt(T(1.0) + r * r)
end
