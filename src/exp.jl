# exported exponential functions

"""
    ldexp(a, n::Int) -> IEEEFloat

Computes `a × 2^n`
"""
ldexp(x::IEEEFloat, q::Int) = ldexpk(x, q)



const over_e2(::Type{Float64}) = 1024
const over_e2(::Type{Float32}) = 128f0

"""
    exp2(x)

Compute the base-`2` exponential of `x`, that is `2ˣ`.
"""
function exp2(x::T) where {T<:IEEEFloat}
    u = expk(dmul(MDLN2(T), x))
    x > over_e2(T) && (u = T(Inf))
    isninf(x) && (u = T(0.0))
    return u
end



const over_e10(::Type{Float64}) = 3.08254715559916743851e2 # log 2^1023*(2-2^-52)
const over_e10(::Type{Float32}) = 38.531839419103626f0 # log 2^127 *(2-2^-23) 

"""
    exp10(x)

Compute the base-`10` exponential of `x`, that is `10ˣ`.
"""
function exp10(x::T) where {T<:IEEEFloat}
    u = expk(dmul(MDLN10(T), x))
    x > over_e10(T) && (u = T(Inf))
    isninf(x) && (u = T(0.0))
    return u
end



const over_em1(::Type{Float64}) = 7.09782712893383996732e2 # log 2^1023*(2-2^-52)
const over_em1(::Type{Float32}) = 88.72283905206835f0 # log 2^127 *(2-2^-23)

const under_em1(::Type{Float64}) = -36.736800569677101399113302437
const under_em1(::Type{Float32}) = -16.635532333438687426013570f0

"""
    expm1(x)

Compute `eˣ- 1` accurately for small values of `x`.
"""
function expm1(x::T) where {T<:IEEEFloat}
    u = T(dadd2(expk2(Double(x)), -T(1.0)))
    x > over_em1(T) && (u = T(Inf))
    x < under_em1(T) && (u = -T(1.0))
    isnegzero(x) && (u = T(-0.0))
    return u
end

const over_e(::Type{Float64}) = 709.78271114955742909217217426
const over_e(::Type{Float32}) = 104f0

const under_e(::Type{Float64}) = -1000.0
const under_e(::Type{Float32}) = -104f0

"""
    exp(x)

Compute the base-`e` exponential of `x`, that is `eˣ`.
"""
function exp end

let
global exp

const c11d = 2.08860621107283687536341e-09
const c10d = 2.51112930892876518610661e-08
const c9d  = 2.75573911234900471893338e-07
const c8d  = 2.75572362911928827629423e-06
const c7d  = 2.4801587159235472998791e-05
const c6d  = 0.000198412698960509205564975
const c5d  = 0.00138888888889774492207962
const c4d  = 0.00833333333331652721664984
const c3d  = 0.0416666666666665047591422
const c2d  = 0.166666666666666851703837
const c1d  = 0.50

const c6f = 0.000198527617612853646278381f0
const c5f = 0.00139304355252534151077271f0
const c4f = 0.00833336077630519866943359f0
const c3f = 0.0416664853692054748535156f0
const c2f = 0.166666671633720397949219f0
const c1f = 0.5f0

global @inline exp_kernel(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d  
global @inline exp_kernel(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f

function exp(d::T) where {T<:IEEEFloat}
    q = unsafe_trunc(Int, round(T(MLN2E) * d))
    s = muladd(q, -L2U(T), d)
    s = muladd(q, -L2L(T), s)

    u = exp_kernel(s)

    u = s * s * u + s + 1
    u = ldexp2k(u, q)

    d < under_e(T) && (u = T(0))
    d > over_e(T) && (u = T(Inf))

    return u
end
end
