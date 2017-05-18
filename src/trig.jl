# exported trigonometric functions

"""
    sin(x)

Compute the sine of `x`, where the output is in radians.
"""
function sin end

"""
    cos(x)

Compute the cosine of `x`, where the output is in radians.
"""
function cos end

let
global sin
global cos

const c8d =  2.72052416138529567917983e-15
const c7d = -7.64292594113954471900203e-13
const c6d =  1.60589370117277896211623e-10
const c5d = -2.5052106814843123359368e-08
const c4d =  2.75573192104428224777379e-06
const c3d = -0.000198412698412046454654947
const c2d =  0.00833333333333318056201922
const c1d = -0.166666666666666657414808

const c4f =  2.6083159809786593541503f-06
const c3f = -0.0001981069071916863322258f0
const c2f =  0.00833307858556509017944336f0
const c1f = -0.166666597127914428710938f0

global @inline _sincos(x::Double{Float64}) = dadd(c1d, x.hi*(@horner x.hi c2d c3d c4d c5d c6d c7d c8d))
global @inline _sincos(x::Double{Float32}) = dadd(c1f, x.hi*(@horner x.hi c2f c3f c4f))

function sin{T<:IEEEFloat}(x::T)
    d = abs(x)
    q = round(d*T(M1PI))
    s = dsub2(d, q*PI4A(T)*4)
    s = dsub2(s, q*PI4B(T)*4)
    s = dsub2(s, q*PI4C(T)*4)
    s = dsub2(s, q*PI4D(T)*4)
    t = s
    s = dsqu(s)
    w =_sincos(s)
    v = dmul(t, dadd(T(1), dmul(w,s)))
    u = T(v)
    qi = unsafe_trunc(Int,q)
    qi & 1 != 0 && (u = -u)
    return flipsign(u,x)
end

function cos{T<:IEEEFloat}(x::T)
    x = abs(x)
    q = muladd(T(2), round(x*T(M1PI) - T(0.5)), T(1))
    s = dsub2(x, q*PI4A(T)*2)
    s = dsub2(s, q*PI4B(T)*2)
    s = dsub2(s, q*PI4C(T)*2)
    s = dsub2(s, q*PI4D(T)*2)
    t = s
    s = dsqu(s)
    w =_sincos(s)
    v = dmul(t, dadd(T(1), dmul(w,s)))
    u = T(v)
    qi = unsafe_trunc(Int,q)
    qi & 2 == 0 && (u = -u)
    return u
end
end


"""
    sin_fast(x)

Compute the sine of `x`, where the output is in radians.
"""
function sin_fast end

"""
    cos_fast(x)

Compute the cosine of `x`, where the output is in radians.
"""
function cos_fast end

let 
global sin_fast
global cos_fast

const c9d = -7.97255955009037868891952e-18
const c8d =  2.81009972710863200091251e-15
const c7d = -7.64712219118158833288484e-13
const c6d =  1.60590430605664501629054e-10
const c5d = -2.50521083763502045810755e-08
const c4d =  2.75573192239198747630416e-06
const c3d = -0.000198412698412696162806809
const c2d =  0.00833333333333332974823815
const c1d = -0.166666666666666657414808

# c5f is 0f0 to handle Inf32 case, Float64 doesn't need this since it comes
# out properly (add another neg constant and remove this zero constant)
const c5f =  0f0
const c4f =  2.6083159809786593541503f-06
const c3f = -0.0001981069071916863322258f0
const c2f =  0.00833307858556509017944336f0
const c1f = -0.166666597127914428710938f0

# Argument is first reduced to the domain 0 < s < π/4

# We return the correct sign using `q & 1 != 0` i.e. q is odd (this works for
# positive and negative q) and if this condition is true we flip the sign since
# we are now in the negative branch of sin(x). Recall that q is just the integer
# part of d/π and thus we can determine the correct sign using this information.

global @inline _sincos_fast(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d
global @inline _sincos_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f

function sin_fast{T<:IEEEFloat}(x::T)
    d = abs(x)
    q = roundi(d*T(M1PI))
    d = muladd(q, -PI4A(T)*4, d)
    d = muladd(q, -PI4B(T)*4, d)
    d = muladd(q, -PI4C(T)*4, d)
    d = muladd(q, -PI4D(T)*4, d)
    s = d*d
    q & 1 != 0 && (d = -d)
    u =_sincos_fast(s)
    u = muladd(s, u*d, d)
    flipsign(u,x)
end

function cos_fast{T<:IEEEFloat}(x::T)
    q = muladd(2, roundi(x*T(M1PI)-T(0.5)), 1)
    x = muladd(q, -PI4A(T)*2, x)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    s = x*x
    q & 2 == 0 && (x = -x)
    u =_sincos_fast(s)
    muladd(s, u*x, x)
end
end



"""
    sincos(x)

Compute the sin and cosine of `x` simultaneously, where the output is in
radians, returning a tuple.
"""
function sincos end

"""
    sincos_fast(x)

Compute the sin and cosine of `x` simultaneously, where the output is in
radians, returning a tuple.
"""
function sincos_fast end

let
global sincos
global sincos_fast

const a6d =  1.58938307283228937328511e-10
const a5d = -2.50506943502539773349318e-08
const a4d =  2.75573131776846360512547e-06
const a3d = -0.000198412698278911770864914
const a2d =  0.0083333333333191845961746
const a1d = -0.166666666666666130709393

const a3f = -0.000195169282960705459117889f0
const a2f =  0.00833215750753879547119141f0
const a1f = -0.166666537523269653320312f0

const b7d = -1.13615350239097429531523e-11
const b6d =  2.08757471207040055479366e-09
const b5d = -2.75573144028847567498567e-07
const b4d =  2.48015872890001867311915e-05
const b3d = -0.00138888888888714019282329
const b2d =  0.0416666666666665519592062
const b1d = -0.50

const b5f = -2.71811842367242206819355f-07
const b4f =  2.47990446951007470488548f-05
const b3f = -0.00138888787478208541870117f0
const b2f =  0.0416666641831398010253906f0
const b1f = -0.5f0

global @inline _sincos_a(x::Float64) = @horner x a1d a2d a3d a4d a5d a6d
global @inline _sincos_a(x::Float32) = @horner x a1f a2f a3f
global @inline _sincos_b(x::Float64) = @horner x b1d b2d b3d b4d b5d b6d b7d
global @inline _sincos_b(x::Float32) = @horner x b1f b2f b3f b4f b5f

function sincos_fast{T<:IEEEFloat}(x::T)
    d  = abs(x)
    q  = roundi(d*T(M2PI))
    s  = d
    s  = muladd(q, -PI4A(T)*2, s)
    s  = muladd(q, -PI4B(T)*2, s)
    s  = muladd(q, -PI4C(T)*2, s)
    s  = muladd(q, -PI4D(T)*2, s)
    t  = s
    s  = s*s
    u  =_sincos_a(s)
    u  = u * s * t
    rx = t + u
    u =_sincos_b(s)
    ry = u * s + T(1)
    q & 1 != 0 && (s = ry; ry = rx; rx = s)
    q & 2 != 0 && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d) && (rx = ry = T(NaN))
    convert(Tuple{T,T}, Double(flipsign(rx,x), ry))
end

function sincos{T<:IEEEFloat}(x::T)
    d  = abs(x)
    q  = roundi(d*2*T(M1PI))
    s  = dsub2(d, q*PI4A(T)*2)
    s  = dsub2(s, q*PI4B(T)*2)
    s  = dsub2(s, q*PI4C(T)*2)
    s  = dsub2(s, q*PI4D(T)*2)
    t  = s
    s  = dsqu(s)
    sx = T(s)
    u  =_sincos_a(sx)
    u *= sx * t.hi
    v  = dadd(t, u)
    rx = T(v)
    u  =_sincos_b(sx)
    v  = dadd(T(1), dmul(sx, u))
    ry = T(v)
    q & 1 != 0 && (u = ry; ry = rx; rx = u)
    q & 2 != 0 && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d) && (rx = ry = T(NaN))
    convert(Tuple{T,T}, Double(flipsign(rx, x), ry))
end
end


"""
    tan(x)

Compute the tangent of `x`, where the output is in radians.
"""
function tan end

"""
    tan_fast(x)

Compute the tangent of `x`, where the output is in radians.
"""
function tan_fast end

let
global tan
global tan_fast

const c15d =  1.01419718511083373224408e-05
const c14d = -2.59519791585924697698614e-05
const c13d =  5.23388081915899855325186e-05
const c12d = -3.05033014433946488225616e-05
const c11d =  7.14707504084242744267497e-05
const c10d =  8.09674518280159187045078e-05
const c9d  =  0.000244884931879331847054404
const c8d  =  0.000588505168743587154904506
const c7d  =  0.001456127889228124279788480
const c6d  =  0.003592087438369066191429240
const c5d  =  0.008863239443624016181133560
const c4d  =  0.021869488285384638959207800
const c3d  =  0.053968253978129841763600200
const c2d  =  0.133333333333125941821962000
const c1d  =  0.333333333333334980164153000

const c7f =  0.00446636462584137916564941f0
const c6f = -8.3920182078145444393158f-05
const c5f =  0.0109639242291450500488281f0
const c4f =  0.0212360303848981857299805f0
const c3f =  0.0540687143802642822265625f0
const c2f =  0.133325666189193725585938f0
const c1f =  0.33333361148834228515625f0

global @inline _tan_fast(x::Float64) = @horner_split x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d
global @inline _tan_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f

function tan_fast{T<:IEEEFloat}(d::T)
    w = abs(d)
    q = roundi(w*T(M2PI))
    x = muladd(q, -PI4A(T)*2, w)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    q & 1 != 0 && (x = -x)
    s = x*x
    u =_tan_fast(s)
    u = muladd(s, u*x, x)
    q & 1 != 0 && (u = 1/u)
    isinf(w) && (u = T(NaN))
    return flipsign(u,d)
end

global @inline _tan(x::Double{Float64}) = dadd(c1d, x.hi*(@horner_split x.hi c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d))
global @inline _tan(x::Double{Float32}) = dadd(c1f, dmul(x, @horner x.hi c2f c3f c4f c5f c6f c7f))
# global @inline _tan(x::Double{Float32}) = dadd(c1f, dmul(x.hi, dadd(c2f, x.hi*(@horner x.hi c3f c4f c5f c6f c7f))))

function tan{T<:IEEEFloat}(d::T)
    w = abs(d)
    q = round(w*T(M2PI))
    x = dsub2(w, q*PI4A(T)*2)
    x = dsub2(x, q*PI4B(T)*2)
    x = dsub2(x, q*PI4C(T)*2)
    x = dsub2(x, q*PI4D(T)*2)
    qi = unsafe_trunc(Int,q)
    qi & 1 != 0 && (x = -x)
    s = dsqu(x)
    u =_tan(s)
    u = dmul(x, dadd(T(1), dmul(u, s)))
    qi & 1 != 0 && (u = ddrec(u))
    return flipsign(T(u),d)
end
end



"""
    atan(x)

Compute the inverse tangent of `x`, where the output is in radians.
"""
function atan{T<:IEEEFloat}(x::T)
    u = T(atan2k(Double(abs(x)), Double(T(1))))
    isinf(x) && (u = T(MPI2))
    flipsign(u,x)
end


"""
    atan_fast(x)

Compute the inverse tangent of `x`, where the output is in radians.
"""
function atan_fast end
let 
global atan_fast

const c19d = -1.88796008463073496563746e-05
const c18d =  0.000209850076645816976906797
const c17d = -0.00110611831486672482563471
const c16d =  0.00370026744188713119232403
const c15d = -0.00889896195887655491740809
const c14d =  0.016599329773529201970117
const c13d = -0.0254517624932312641616861
const c12d =  0.0337852580001353069993897
const c11d = -0.0407629191276836500001934
const c10d =  0.0466667150077840625632675
const c9d  = -0.0523674852303482457616113
const c8d  =  0.0587666392926673580854313
const c7d  = -0.0666573579361080525984562
const c6d  =  0.0769219538311769618355029
const c5d  = -0.090908995008245008229153
const c4d  =  0.111111105648261418443745
const c3d  = -0.14285714266771329383765
const c2d  =  0.199999999996591265594148
const c1d  = -0.333333333333311110369124

const c8f =  0.00282363896258175373077393f0
const c7f = -0.0159569028764963150024414f0
const c6f =  0.0425049886107444763183594f0
const c5f = -0.0748900920152664184570312f0
const c4f =  0.106347933411598205566406f0
const c3f = -0.142027363181114196777344f0
const c2f =  0.199926957488059997558594f0
const c1f = -0.333331018686294555664062f0

global @inline _atan_fast(x::Float64) = @horner_split x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d
global @inline _atan_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f

function atan_fast{T<:IEEEFloat}(x::T)
    q = 0
    if x < 0
        x = -x
        q = 2
    end
    if x > 1
        x = 1/x
        q |= 1
    end
    t = x*x
    u =_atan_fast(t)
    t = x + x*t*u
    q & 1 != 0 && (t = T(MPI2) - t)
    q & 2 != 0 && (t = -t)
    return t
end
end



"""
    atan2(x, y)

Compute the inverse tangent of `x/y`, using the signs of both `x` and `y` to determine the quadrant of the return value.
"""
function atan2{T<:IEEEFloat}(x::T, y::T)
    r = T(atan2k(Double(abs(x)), Double(y)))
    r = flipsign(r,y)
    if isinf(y) || y == 0
        r = T(MPI2) - (isinf(y) ? _sign(y)*T(MPI2) : T(0))
    end
    if isinf(x)
        r = T(MPI2) - (isinf(y) ? _sign(y)*T(MPI4) : T(0))
    end
    if x == 0
        r = _sign(y) == -1 ? T(MPI) : T(0)
    end
    return isnan(y) || isnan(x) ? T(NaN) : flipsign(r,x)
end


"""
    atan2_fast(x, y)

Compute the inverse tangent of `x/y`, using the signs of both `x` and `y` to determine the quadrant of the return value.
"""
function atan2_fast{T<:IEEEFloat}(x::T, y::T)
    r = atan2k_fast(abs(x), y)
    r = flipsign(r,y)
    if isinf(y) || y == 0
        r = T(MPI2) - (isinf(y) ? _sign(y)*T(MPI2) : T(0))
    end
    if isinf(x)
        r = T(MPI2) - (isinf(y) ? _sign(y)*T(MPI4) : T(0))
    end
    if x == 0
        r = _sign(y) == -1 ? T(MPI) : T(0)
    end
    return isnan(y) || isnan(x) ? T(NaN) : flipsign(r,x)
end



"""
    asin(x)

Compute the inverse sine of `x`, where the output is in radians.
"""
function asin{T<:IEEEFloat}(x::T)
    d = atan2k(Double(abs(x)), dsqrt(dmul(dadd(T(1), x), dsub(T(1), x))))
    u = T(d)
    abs(x) == 1 && (u = T(MPI2))
    flipsign(u,x)
end


"""
    asin_fast(x)

Compute the inverse sine of `x`, where the output is in radians.
"""
function asin_fast{T<:IEEEFloat}(x::T)
    flipsign(atan2k_fast(abs(x), _sqrt((1+x)*(1-x))), x)
end



"""
    acos(x)

Compute the inverse cosine of `x`, where the output is in radians.
"""
function acos{T<:IEEEFloat}(x::T)
    d = atan2k(dsqrt(dmul(dadd(T(1), x), dsub(T(1), x))), Double(abs(x)))
    d = flipsign(d,x)
    abs(x) == 1 && (d = Double(T(0)))
    x < 0 && (d = dadd(MDPI(T), d))
    return T(d)
end


"""
    acos_fast(x)

Compute the inverse cosine of `x`, where the output is in radians.
"""
function acos_fast{T<:IEEEFloat}(x::T)
    flipsign(atan2k_fast(_sqrt((1+x)*(1-x)), abs(x)), x) + (x < 0 ? T(MPI) : T(0))
end
