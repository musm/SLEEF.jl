# exported exponential functions

"""
    ldexp(a, n::Int) -> IEEEFloat

Computes `a × 2^n`
"""
ldexp(x::IEEEFloat, q::Int) = ldexpk(x,q)



const over_e2(::Type{Float64}) = 1024
const over_e2(::Type{Float32}) = 128f0

"""
    exp2(x)

Compute the base-`2` exponential of `x`, that is `2ˣ`.
"""
function exp2{T<:IEEEFloat}(x::T)
    u = expk(dmul(MDLN2(T), x))
    x > over_e2(T) && (u = T(Inf))
    isninf(x) && (u = T(0))
    return u
end



const over_e10(::Type{Float64}) = 308
const over_e10(::Type{Float32}) = 38f0

"""
    exp10(x)

Compute the base-`10` exponential of `x`, that is `10ˣ`.
"""
function exp10{T<:IEEEFloat}(x::T)
    u = expk(dmul(MDLN10(T), x))
    x > over_e10(T) && (u = T(Inf))
    isninf(x) && (u = T(0))
    return u
end



const over_em1(::Type{Float64}) = 700.0
const over_em1(::Type{Float32}) = 88f0

const under_em1(::Type{Float64}) = -0.36043653389117156089696070315825181539851971360337e2
const under_em1(::Type{Float32}) = -0.15942385152878742116596338793538061065739925620174f2

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

const c5f = 0.00136324646882712841033936f0
const c4f = 0.00836596917361021041870117f0
const c3f = 0.0416710823774337768554688f0
const c2f = 0.166665524244308471679688f0
const c1f = 0.499999850988388061523438f0

global @inline _exp(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d  
global @inline _exp(x::Float32) = @horner x c1f c2f c3f c4f c5f

function exp{T<:IEEEFloat}(x::T)
    q = roundi(T(MLN2E)*x)
    s = muladd(q, -LN2U(T), x)
    s = muladd(q, -LN2L(T), s)
    u =_exp(s)
    u = s*s*u + s + 1
    u = ldexpk(u,q)
    isninf(x) && (u = T(0))
    return u
end
end
