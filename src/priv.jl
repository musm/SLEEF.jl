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
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    m, q
end
@inline split_exponent(::Type{Float64}, q::Int) = _split_exponent(q, UInt(9), UInt(31), UInt(2))
@inline split_exponent(::Type{Float32}, q::Int) = _split_exponent(q, UInt(6), UInt(31), UInt(2))

"""
    ldexpk(a::IEEEFloat, n::Int) -> IEEEFloat

Computes `a × 2^n`.
"""
@inline function ldexpk{T<:IEEEFloat}(x::T, q::Int)
    bias = exponent_bias(T)
    emax = exponent_max(T)
    m, q = split_exponent(T,q)
    m += bias
    m = ifelse(m < 0, 0, m)
    m = ifelse(m > emax, emax, m)
    q += bias
    u = integer2float(T, m)
    x = x*u*u*u*u
    u = integer2float(T, q)
    x*u
end

@inline ldexpk2{T<:IEEEFloat}(x::T, e::Int) = x * pow2i(T, e >> UInt(1)) * pow2i(T, e - (e >> UInt(1)))



# The following define threshold values for `ilog2k`
threshold_exponent(::Type{Float64}) = 300
threshold_exponent(::Type{Float32}) = 64

"""
    ilog2k(x::IEEEFloat) -> Int

Returns the integral part of the logarithm of `|x|`, using 2 as base for the logarithm; in other
words this returns the binary exponent of `x` so that

    x = significand × 2^exponent

where `significand ∈ [0.5, 1)`.
"""
@inline function ilog2k{T<:IEEEFloat}(d::T)
    m = d < T(2)^-threshold_exponent(T)
    d = ifelse(m, d*T(2)^threshold_exponent(T), d)
    q = float2integer(d) & exponent_max(T)
    q = ifelse(m, q - (threshold_exponent(T) + exponent_bias(T) - 1), q - (exponent_bias(T) - 1)) # subtract 1 since we want 2^q
end


let
const c20d =  1.06298484191448746607415e-05
const c19d = -0.000125620649967286867384336
const c18d =  0.00070557664296393412389774
const c17d = -0.00251865614498713360352999
const c16d =  0.00646262899036991172313504
const c15d = -0.0128281333663399031014274
const c14d =  0.0208024799924145797902497
const c13d = -0.0289002344784740315686289
const c12d =  0.0359785005035104590853656
const c11d = -0.041848579703592507506027
const c10d =  0.0470843011653283988193763
const c9d  = -0.0524914210588448421068719
const c8d  =  0.0587946590969581003860434
const c7d  = -0.0666620884778795497194182
const c6d  =  0.0769225330296203768654095
const c5d  = -0.0909090442773387574781907
const c4d  =  0.111111108376896236538123
const c3d  = -0.142857142756268568062339
const c2d  =  0.199999999997977351284817
const c1d =  -0.333333333333317605173818

const c9f = -0.00176397908944636583328247f0
const c8f =  0.0107900900766253471374512f0
const c7f = -0.0309564601629972457885742f0
const c6f =  0.0577365085482597351074219f0
const c5f = -0.0838950723409652709960938f0
const c4f =  0.109463557600975036621094f0
const c3f = -0.142626821994781494140625f0
const c2f =  0.199983194470405578613281f0
const c1f = -0.333332866430282592773438f0

global @inline _atan2k_fast(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline _atan2k_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f c9f

global @inline function atan2k_fast{T<:IEEEFloat}(y::T, x::T)
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
    s = y/x
    t = s*s
    u =_atan2k_fast(t)
    t = u*t*s + s
    q*T(MPI2) + t
end


global @inline _atan2k(x::Double{Float64}) = @horner x.hi c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline _atan2k(x::Double{Float32}) = dadd(c1f, x.hi*(@horner x.hi c2f c3f c4f c5f c6f c7f c8f c9f))

global @inline function atan2k{T<:IEEEFloat}(y::Double{T}, x::Double{T})
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
    s = ddiv(y,x)
    t = dsqu(s)
    u =_atan2k(t)
    t = dmul(s, dadd(T(1), dmul(t, u)))
    dadd(dmul(T(q), MDPI2(T)), t)
end
end


let
const c10d = 2.51069683420950419527139e-08
const c9d  = 2.76286166770270649116855e-07
const c8d  = 2.75572496725023574143864e-06
const c7d  = 2.48014973989819794114153e-05
const c6d  = 0.000198412698809069797676111
const c5d  = 0.0013888888939977128960529
const c4d  = 0.00833333333332371417601081
const c3d  = 0.0416666666665409524128449
const c2d  = 0.166666666666666740681535
const c1d  = 0.500000000000000999200722

const c5f = 0.00136324646882712841033936f0
const c4f = 0.00836596917361021041870117f0
const c3f = 0.0416710823774337768554688f0
const c2f = 0.166665524244308471679688f0
const c1f = 0.499999850988388061523438f0

global @inline _expk(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d
global @inline _expk(x::Float32) = @horner x c1f c2f c3f c4f c5f

global @inline function expk{T<:IEEEFloat}(d::Double{T})
    q = round(T(d)*T(MLN2E))
    s = dadd(d, -q*LN2U(T))
    s = dadd(s, -q*LN2L(T))
    u =_expk(T(s))
    t = dadd(s, dmul(dsqu(s), u))
    t = dadd(T(1), t)
    ldexpk(T(t), unsafe_trunc(Int,q))
end


global @inline function expk2{T<:IEEEFloat}(d::Double{T})
    q = round(T(d)*T(MLN2E))
    s = dadd(d, -q*LN2U(T))
    s = dadd(s, -q*LN2L(T))
    u =_expk(s.hi)
    t = dadd(s, dmul(dsqu(s), u))
    t = dadd(T(1), t)
    scale(t, pow2i(T, unsafe_trunc(Int,q)))
end
end


let
const c8d = 0.134601987501262130076155
const c7d = 0.132248509032032670243288
const c6d = 0.153883458318096079652524
const c5d = 0.181817427573705403298686
const c4d = 0.222222231326187414840781
const c3d = 0.285714285651261412873718
const c2d = 0.400000000000222439910458
const c1d = 0.666666666666666371239645

const c4f = 0.2371599674224853515625f0
const c3f = 0.285279005765914916992188f0
const c2f = 0.400005519390106201171875f0
const c1f = 0.666666567325592041015625f0

global @inline _logk(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d
global @inline _logk(x::Float32) = @horner x c1f c2f c3f c4f

global @inline function logk{T<:IEEEFloat}(d::T)
    e  = ilog2k(d*T(M1SQRT2))
    m  = ldexpk(d,-e)
    x  = ddiv(dsub2(m, T(1)), dadd2(T(1), m))
    x2 = dsqu(x)
    t  =_logk(x2.hi)
    dadd(dmul(MDLN2(T), T(e)), dadd(scale(x, T(2)), dmul(dmul(x2, x), t)))
end


global @inline function logk2{T<:IEEEFloat}(d::Double{T})
    e  = ilog2k(d.hi*T(M1SQRT2))
    m  = scale(d, pow2i(T,-e))
    x  = ddiv(dsub2(m, T(1)), dadd2(m, T(1)))
    x2 = dsqu(x)
    t  =_logk(x2.hi)
    dadd(dmul(MDLN2(T), T(e)), dadd(scale(x, T(2)), dmul(dmul(x2, x), t)))
end
end
