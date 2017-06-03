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
function log10(a::T) where {T<:IEEEFloat}
    x = T(dmul(logk(a), MDLN10E(T)))

    isinf(a) && (x = T(Inf))
    a < 0  && (x = T(NaN))
    a == 0 && (x = T(-Inf))

    return x
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


const over_log1p(::Type{Float64}) = 1e307
const over_log1p(::Type{Float32}) = 1f38

"""
    log1p(x)

Accurately compute the natural logarithm of 1+x.
"""
function log1p(a::T) where {T<:IEEEFloat}
    x = T(logk2(dadd2(a, T(1.0))))

    a > over_log1p(T) && (x = T(Inf))
    a < -1 && (x = T(NaN))
    a == -1 && (x = T(-Inf))
    isnegzero(a) && (x = T(-0.0))

    return x
end


"""
    log(x)

Compute the natural logarithm of `x`. The inverse of the natural logarithm is
the natural expoenential function `exp(x)`
"""
function log(d::T) where {T<:IEEEFloat}
    x = T(logk(d))

    isinf(d) && (x = T(Inf))
    d < 0  && (x = T(NaN))
    d == 0 && (x = -T(Inf))

    return x
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

const c8d = 0.153487338491425068243146
const c7d = 0.152519917006351951593857
const c6d = 0.181863266251982985677316
const c5d = 0.222221366518767365905163
const c4d = 0.285714294746548025383248
const c3d = 0.399999999950799600689777
const c2d = 0.6666666666667778740063
const c1d = 2.0

const c5f = 0.2392828464508056640625f0
const c4f = 0.28518211841583251953125f0
const c3f = 0.400005877017974853515625f0
const c2f = 0.666666686534881591796875f0
const c1f = 2f0

global @inline log_fast_kernel(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d
global @inline log_fast_kernel(x::Float32) = @horner x c1f c2f c3f c4f c5f

function log_fast(d::T) where {T<:IEEEFloat}
    e  = ilogbk(d * T(1.0/0.75))
    m  = ldexpk(d, -e)

    x  = (m - 1) / (m + 1)
    x2 = x * x

    t  = log_fast_kernel(x2)
    
    x  = x * t + T(MLN2) * e
    
    isinf(d) && (x = T(Inf))
    d < 0 && (x = T(NaN))
    d == 0 && (x = -T(Inf))

    return x
end
end
