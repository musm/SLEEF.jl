# exported logarithmic functions
const FP_ILOGB0   = typemin(Int)
const FP_ILOGBNAN = typemin(Int)
const INT_MAX     = typemax(Int)
"""
    ilogb(x)

Returns the integral part of the logarithm of `abs(x)`, using base 2 for the
logarithm. In other words, this computes the binary exponent of `x` such that

    x = significand × 2^exponent,

where `significand ∈ [1, 2)`.

* Exceptional cases (where `Int` is the machine wordsize)
    * `x = 0`    returns `FP_ILOGB0`
    * `x = ±Inf`  returns `INT_MAX`
    * `x = NaN`  returns `FP_ILOGBNAN`
"""
function ilogb(x::T) where {T<:IEEEFloat}
    e = ilogbk(abs(x))
    x == 0 && (e = FP_ILOGB0)
    isnan(x) && (e = FP_ILOGBNAN)
    isinf(x) && (e = INT_MAX)
    return e
end


"""
    log10(x)

Returns the base `10` logarithm of `x`.
"""
function log10(x::T) where {T<:IEEEFloat}
    u = T(dmul(logk(x), MDLN10E(T)))
    isinf(x) && (u = T(Inf))
    x < 0  && (u = T(NaN))
    x == 0 && (u = T(-Inf))
    return u
end


"""
    log2(x)

Returns the base `2` logarithm of `x`.
"""
function log2(x::T) where {T<:IEEEFloat}
    u = T(dmul(logk(x), MDLN2E(T)))
    isinf(x) && (u = T(Inf))
    x < 0  && (u = T(NaN))
    x == 0 && (u = T(-Inf))
    return u
end


"""
    log1p(x)

Accurately compute the natural logarithm of 1+x.
"""
function log1p(x::T) where {T<:IEEEFloat}
    u = T(logk2(dadd2(x, T(1.0))))

    isinf(x) && (u = T(Inf))
    x < -1 && (u = T(NaN))
    x == -1 && (u = -T(Inf))
    isnegzero(x) && (u = T(-0.0))

    return u
end


"""
    log(x)

Compute the natural logarithm of `x`. The inverse of the natural logarithm is
the natural expoenential function `exp(x)`
"""
function log(x::T) where {T<:IEEEFloat}
    u = T(logk(x))
    isinf(x) && (u = T(Inf))
    x < 0  && (u = T(NaN))
    x == 0 && (u = -T(Inf))
    return u
end

# First we split the argument to its mantissa `m` and integer exponent `e` so
# that `d = m \times 2^e`, where `m \in [0.5, 1)` then we apply the polynomial
# approximant on this reduced argument `m` before putting back the exponent
# in. This first part is done with the help of the private function
# `ilogbk(x)` and we put the exponent back using

#     `\log(m \times 2^e) = \log(m) + \log 2^e =  \log(m) + e\times MLN2

# The polynomial we evaluate is based on coefficients from

#     `log_2(x) = 2\sum_{n=0}^\infty \frac{1}{2n+1} \bigl(\frac{x-1}{x+1}^{2n+1}\bigr)`

# That being said, since this converges faster when the argument is close to
# 1, we multiply  `m` by `2` and subtract 1 for the exponent `e` when `m` is
# less than `sqrt(2)/2`
"""
    log_fast(x)

Compute the natural logarithm of `x`. The inverse of the natural logarithm is
the natural expoenential function `exp(x)`
"""
function log_fast end

let
global log_fast

const c8d = 0.148197055177935105296783
const c7d = 0.153108178020442575739679
const c6d = 0.181837339521549679055568
const c5d = 0.22222194152736701733275
const c4d = 0.285714288030134544449368
const c3d = 0.399999999989941956712869
const c2d = 0.666666666666685503450651
const c1d = 2.0

const c5f = 0.2371599674224853515625f0
const c4f = 0.285279005765914916992188f0
const c3f = 0.400005519390106201171875f0
const c2f = 0.666666567325592041015625f0
const c1f = 2f0

global @inline log_fast_kernel(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d
global @inline log_fast_kernel(x::Float32) = @horner x c1f c2f c3f c4f c5f

function log_fast(x::T) where {T<:IEEEFloat}
    e  = ilogbk(T(MSQRT2)*x)
    m  = ldexpk(x, -e)
    u  = (m-1) / (m+1)
    u2 = u*u

    t  = log_fast_kernel(u2)
    u  = muladd(u, t, T(MLN2)*e)
    isinf(x) && (u = T(Inf))
    x < 0 && (u = T(NaN))
    x == 0 && (u = -T(Inf))
    return u
end
end
