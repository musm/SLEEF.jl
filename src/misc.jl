
"""
  pow(x, y)

Exponentiation operator, returns `x` raised to the power `y`.
"""
function pow(x::V, y::V) where {V <: FloatType}
    T = eltype(x)
    yi = unsafe_trunc(fpinttype(T), y)
    yisint = yi == y
    yisodd = isodd(yi) & yisint

    result = expk(dmul(logk(abs(x)), y))

    result = vifelse(isnan(result), T(Inf), result)
    result = vifelse(x > 0, result, vifelse(!yisint, T(NaN), vifelse(yisodd, -result, result)))

    efx = flipsign(abs(x) - 1, y)
    result = vifelse(isinf(y), vifelse(efx < 0, T(0.0), vifelse(efx == 0, T(1.0), T(Inf))), result)
    result = vifelse(isinf(x) | (x == 0), vifelse(yisodd, _sign(x), T(1.0)) * vifelse(vifelse(x == 0, -y, y) < 0, T(0.0), T(Inf)), result)
    result = vifelse(isnan(x) | isnan(y), T(NaN), result)
    result = vifelse((y == 0) | (x == 1), T(1.0), result)

    return result
end






@inline function cbrt_kernel(x::FloatType64)
    c6d = -0.640245898480692909870982
    c5d = 2.96155103020039511818595
    c4d = -5.73353060922947843636166
    c3d = 6.03990368989458747961407
    c2d = -3.85841935510444988821632
    c1d = 2.2307275302496609725722
    @horner x c1d c2d c3d c4d c5d c6d
end
@inline function cbrt_kernel(x::FloatType32)
    c6f = -0.601564466953277587890625f0
    c5f =  2.8208892345428466796875f0
    c4f = -5.532182216644287109375f0
    c3f =  5.898262500762939453125f0
    c2f = -3.8095417022705078125f0
    c1f =  2.2241256237030029296875f0
    @horner x c1f c2f c3f c4f c5f c6f
end

"""
    cbrt_fast(x)

Return `x^{1/3}`.
"""
function cbrt_fast(d::V) where {V <: FloatType}
    T  = eltype(d)
    e  = ilogbk(abs(d)) + 1
    d  = ldexp2k(d, -e)
    r  = (e + 6144) % 3
    q  = vifelse(r == 1, V(M2P13), V(1))
    q  = vifelse(r == 2, V(M2P23), q)
    q  = ldexp2k(q, (e + 6144) ÷ 3 - 2048)
    q  = flipsign(q, d)
    d  = abs(d)
    x  = cbrt_kernel(d)
    y  = x * x
    y  = y * y
    x -= (d * y - x) * T(1 / 3)
    y  = d * x * x
    y = (y - T(2 / 3) * y * (y * x - 1)) * q
end


"""
    cbrt(x)

Return `x^{1/3}`. The prefix operator `∛` is equivalent to `cbrt`.
"""
function cbrt(d::V) where {V <: FloatType}
    T  = eltype(d)
    e  = ilogbk(abs(d)) + 1
    d  = ldexp2k(d, -e)
    r  = (e + 6144) % 3
    q2 = vifelse(r == 1, MD2P13(T), Double(T(1)))
    q2 = vifelse(r == 2, MD2P23(T), q2)
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
    y  = V(u)
    y  = T(-2 / 3) * y * z
    v  = dadd(dsqu(z), y)
    v  = dmul(v, d)
    v  = dmul(v, q2)
    z  = ldexp2k(V(v), (e + 6144) ÷ 3 - 2048)
    z  = vifelse(isinf(d), flipsign(T(Inf), q2.hi), z)
    z  = vifelse(d == 0, flipsign(T(0), q2.hi), z)
    return z
end



"""
    hypot(x,y)

Compute the hypotenuse `\\sqrt{x^2+y^2}` avoiding overflow and underflow.
"""
function hypot(x::T, y::T) where {T<:vIEEEFloat}
    a = abs(x)
    b = abs(y)

    x = min(a,b)
    y = max(a,b)

    r = vifelse(x == 0, y, y / x)
    x * sqrt(T(1.0) + r * r)
end
