# private math functions

"""
A helper function for `ldexpk`

First note that `r = (q >> n) << n` clears the lowest n bits of q, i.e. returns 2^n where n is the
largest integer such that q >= 2^n

For numbers q less than 2^m the following code does the same as the above snippet
    `r = ( (q>>v + q) >> n - q>>v ) << n`
For numbers larger than or equal to 2^v this subtracts 2^n from q for q>>n times.

The function returns q(input) := q(output) + offset*r

In the code for ldexpk we actually use
    `m = ( (m>>n + m) >> n -  m>>m ) << (n-2)`.
So that x has to be multplied by u four times `x = x*u*u*u*u` to put the value  of the offset
exponent amount back in.
"""
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n - offset)
    q = q - (m << offset)
    m, q
end
@inline split_exponent(::Type{Float64}, q::Int) = _split_exponent(q, UInt(9), UInt(31), UInt(2))
@inline split_exponent(::Type{Float32}, q::Int) = _split_exponent(q, UInt(6), UInt(31), UInt(2))

"""
    ldexpk(a, n)

Computes `a × 2^n`.
"""
@inline function ldexpk(x::T, q::Int) where {T<:Union{Float32,Float64}}
    bias = exponent_bias(T)
    emax = exponent_raw_max(T)
    m, q = split_exponent(T, q)
    m += bias
    m = ifelse(m < 0, 0, m)
    m = ifelse(m > emax, emax, m)
    q += bias
    u = integer2float(T, m)
    x = x * u * u * u * u
    u = integer2float(T, q)
    x * u
end

@inline function ldexp2k(x::T, e::Int) where {T<:Union{Float32,Float64}}
    x * pow2i(T, e >> 1) * pow2i(T, e - (e >> 1))
end

@inline function ldexp3k(x::T, e::Int) where {T<:Union{Float32,Float64}}
    reinterpret(T, reinterpret(Unsigned, x) + (Int64(e) << significand_bits(T)) % uinttype(T))
end

# threshold values for `ilogbk`
const threshold_exponent(::Type{Float64}) = 300
const threshold_exponent(::Type{Float32}) = 64

"""
    ilogbk(x) -> Int

Returns the integral part of the logarithm of `|x|`, using 2 as base for the logarithm; in other
words this returns the binary exponent of `x` so that

    x = significand × 2^exponent

where `significand ∈ [1, 2)`.
"""
@inline function ilogbk(d::T) where {T<:Union{Float32,Float64}}
    m = d < T(2)^-threshold_exponent(T)
    d = ifelse(m, d * T(2)^threshold_exponent(T), d)
    q = float2integer(d) & exponent_raw_max(T)
    q = ifelse(m, q - (threshold_exponent(T) + exponent_bias(T)), q - exponent_bias(T))
end

# similar to ilogbk, but argument has to be a normalized float value
@inline function ilogb2k(d::T) where {T<:Union{Float32,Float64}}
    (float2integer(d) & exponent_raw_max(T)) - exponent_bias(T)
end


let
global atan2k_fast
global atan2k

c20d =  1.06298484191448746607415e-05
c19d = -0.000125620649967286867384336
c18d =  0.00070557664296393412389774
c17d = -0.00251865614498713360352999
c16d =  0.00646262899036991172313504
c15d = -0.0128281333663399031014274
c14d =  0.0208024799924145797902497
c13d = -0.0289002344784740315686289
c12d =  0.0359785005035104590853656
c11d = -0.041848579703592507506027
c10d =  0.0470843011653283988193763
c9d  = -0.0524914210588448421068719
c8d  =  0.0587946590969581003860434
c7d  = -0.0666620884778795497194182
c6d  =  0.0769225330296203768654095
c5d  = -0.0909090442773387574781907
c4d  =  0.111111108376896236538123
c3d  = -0.142857142756268568062339
c2d  =  0.199999999997977351284817
c1d =  -0.333333333333317605173818

c9f = -0.00176397908944636583328247f0
c8f =  0.0107900900766253471374512f0
c7f = -0.0309564601629972457885742f0
c6f =  0.0577365085482597351074219f0
c5f = -0.0838950723409652709960938f0
c4f =  0.109463557600975036621094f0
c3f = -0.142626821994781494140625f0
c2f =  0.199983194470405578613281f0
c1f = -0.333332866430282592773438f0

global @inline atan2k_fast_kernel(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline atan2k_fast_kernel(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f c9f

@inline function atan2k_fast(y::T, x::T) where {T<:Union{Float32,Float64}}
    q = 0
    if x < 0
        x = -x
        q = -2
    end
    if y > x
        t = x; x = y
        y = -t
        q += 1
    end
    s = y / x
    t = s * s
    u = atan2k_fast_kernel(t)
    t = u * t * s + s
    q * T(PI_2) + t
end


global @inline atan2k_kernel(x::Double{Float64}) = @horner x.hi c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline atan2k_kernel(x::Double{Float32}) = dadd(c1f, x.hi * (@horner x.hi c2f c3f c4f c5f c6f c7f c8f c9f))

@inline function atan2k(y::Double{T}, x::Double{T}) where {T<:Union{Float32,Float64}}
    q = 0
    if x < 0
        x = -x
        q = -2
    end
    if y > x
        t = x; x = y
        y = -t
        q += 1
    end

    s = ddiv(y, x)
    t = dsqu(s)
    t = dnormalize(t)

    u = atan2k_kernel(t)

    t = dmul(t, u)
    t = dmul(s, dadd(T(1.0), t))
    T <: Float64 && abs(s.hi) < 1e-200 && (t = s)
    t = dadd(dmul(T(q), MDPI2(T)), t)
    return t
end
end



const under_expk(::Type{Float64}) = -1000.0
const under_expk(::Type{Float32}) = -104f0

@inline function expk_kernel(x::Float64)
    c10 = 2.51069683420950419527139e-08
    c9  = 2.76286166770270649116855e-07
    c8  = 2.75572496725023574143864e-06
    c7  = 2.48014973989819794114153e-05
    c6  = 0.000198412698809069797676111
    c5  = 0.0013888888939977128960529
    c4  = 0.00833333333332371417601081
    c3  = 0.0416666666665409524128449
    c2  = 0.166666666666666740681535
    c1  = 0.500000000000000999200722
    return @horner x c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
end

@inline function  expk_kernel(x::Float32)
    c5 = 0.00136324646882712841033936f0
    c4 = 0.00836596917361021041870117f0
    c3 = 0.0416710823774337768554688f0
    c2 = 0.166665524244308471679688f0
    c1 = 0.499999850988388061523438f0
    return @horner x c1 c2 c3 c4 c5
end

@inline function expk(d::Double{T}) where {T<:Union{Float32,Float64}}
    q = round(T(d) * T(MLN2E))
    qi = unsafe_trunc(Int, q)

    s = dadd(d, -q * L2U(T))
    s = dadd(s, -q * L2L(T))

    s = dnormalize(s)

    u = expk_kernel(T(s))

    t = dadd(s, dmul(dsqu(s), u))
    t = dadd(T(1.0), t)
    u = ldexpk(T(t), qi)

    (d.hi < under_expk(T)) && (u = T(0.0))
    return u
end



@inline function expk2_kernel(x::Double{Float64})
    c11 = 0.1602472219709932072e-9
    c10  = 0.2092255183563157007e-8
    c9  = 0.2505230023782644465e-7
    c8  = 0.2755724800902135303e-6
    c7  = 0.2755731892386044373e-5
    c6  = 0.2480158735605815065e-4
    c5  = 0.1984126984148071858e-3
    c4  = 0.1388888888886763255e-2
    c3  = 0.8333333333333347095e-2
    c2  = 0.4166666666666669905e-1
    c1 = 0.1666666666666666574e0
    u = @horner x.hi c2 c3 c4 c5 c6 c7 c8 c9 c10 c11
    return dadd(dmul(x, u), c1)
end

@inline function  expk2_kernel(x::Double{Float32})
    c5 = 0.1980960224f-3
    c4 = 0.1394256484f-2
    c3 = 0.8333456703f-2
    c2 = 0.4166637361f-1
    c1 = 0.166666659414234244790680580464f0
    u = @horner x.hi c2 c3 c4 c5
    return dadd(dmul(x, u), c1)
end

@inline function expk2(d::Double{T}) where {T<:Union{Float32,Float64}}
    q = round(T(d) * T(MLN2E))
    qi = unsafe_trunc(Int, q)

    s = dadd(d, -q * L2U(T))
    s = dadd(s, -q * L2L(T))

    t = expk2_kernel(s)
    t = dadd(dmul(s, t), T(0.5))
    t = dadd(s, dmul(dsqu(s), t))
    t = dadd(T(1.0), t)

    t = Double(ldexp2k(t.hi, qi), ldexp2k(t.lo, qi))

    (d.hi < under_expk(T)) && (t = Double(T(0.0)))
    return t
end



@inline function logk2_kernel(x::Float64)
    c8 = 0.13860436390467167910856
    c7 = 0.131699838841615374240845
    c6 = 0.153914168346271945653214
    c5 = 0.181816523941564611721589
    c4 = 0.22222224632662035403996
    c3 = 0.285714285511134091777308
    c2 = 0.400000000000914013309483
    c1 = 0.666666666666664853302393
    return @horner x c1 c2 c3 c4 c5 c6 c7 c8
end

@inline function logk2_kernel(x::Float32)
    c4 = 0.240320354700088500976562f0
    c3 = 0.285112679004669189453125f0
    c2 = 0.400007992982864379882812f0
    c1 = 0.666666686534881591796875f0
    return @horner x c1 c2 c3 c4
end

@inline function logk2(d::Double{T}) where {T<:Union{Float32,Float64}}
    e  = ilogbk(d.hi * T(1.0/0.75))
    m  = scale(d, pow2i(T, -e))

    x  = ddiv(dsub2(m, T(1.0)), dadd2(m, T(1.0)))
    x2 = dsqu(x)

    t  = logk2_kernel(x2.hi)

    s = dmul(MDLN2(T), T(e))
    s = dadd(s, scale(x, T(2.0)))
    s = dadd(s, dmul(dmul(x2, x), t))
    return s
end



@inline function logk_kernel(x::Double{Float64})
    c10 = 0.116255524079935043668677
    c9 = 0.103239680901072952701192
    c8 = 0.117754809412463995466069
    c7 = 0.13332981086846273921509
    c6 = 0.153846227114512262845736
    c5 = 0.181818180850050775676507
    c4 = 0.222222222230083560345903
    c3 = 0.285714285714249172087875
    c2 = 0.400000000000000077715612
    c1 = Double(0.666666666666666629659233, 3.80554962542412056336616e-17)
    dadd2(dmul(x, @horner x.hi c2 c3 c4 c5 c6 c7 c8 c9 c10), c1)
end

@inline function logk_kernel(x::Double{Float32})
    c4 = 0.240320354700088500976562f0
    c3 = 0.285112679004669189453125f0
    c2 = 0.400007992982864379882812f0
    c1 = Double(0.66666662693023681640625f0, 3.69183861259614332084311f-9)    
    dadd2(dmul(x, @horner x.hi c2 c3 c4), c1)
end

@inline function logk(d::T) where {T<:Union{Float32,Float64}}
    o = d < floatmin(T)
    o && (d *= T(Int64(1) << 32) * T(Int64(1) << 32))

    e  = ilogb2k(d * T(1.0/0.75))
    m  = ldexp3k(d, -e)

    o && (e -= 64)

    x  = ddiv(dsub2(m, T(1.0)), dadd2(T(1.0), m))
    x2 = dsqu(x)

    t  = logk_kernel(x2)

    s = dmul(MDLN2(T), T(e))
    s = dadd(s, scale(x, T(2.0)))
    s = dadd(s, dmul(dmul(x2, x), t))
    return s
end
