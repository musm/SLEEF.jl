# exported exponential functions

"""
    ldexp(a, n)

Computes `a × 2^n`
"""
ldexp(x::Union{Float32,Float64}, q::Int) = ldexpk(x, q)


const max_exp2(::Type{Float64}) = 1024
const max_exp2(::Type{Float32}) = 128f0

const min_exp2(::Type{Float64}) = -1075
const min_exp2(::Type{Float32}) = -150f0

@inline function exp2_kernel(x::Float64)
    c11 = 0.4434359082926529454e-9
    c10 = 0.7073164598085707425e-8
    c9  = 0.1017819260921760451e-6
    c8  = 0.1321543872511327615e-5
    c7  = 0.1525273353517584730e-4
    c6  = 0.1540353045101147808e-3
    c5  = 0.1333355814670499073e-2
    c4  = 0.9618129107597600536e-2
    c3  = 0.5550410866482046596e-1
    c2  = 0.2402265069591012214
    c1  = 0.6931471805599452862
    return @horner x c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11
end

@inline function exp2_kernel(x::Float32)
    c6 = 0.1535920892f-3
    c5 = 0.1339262701f-2
    c4 = 0.9618384764f-2
    c3 = 0.5550347269f-1
    c2 = 0.2402264476f0
    c1 = 0.6931471825f0
    return @horner x c1 c2 c3 c4 c5 c6
end

"""
    exp2(x)

Compute the base-`2` exponential of `x`, that is `2ˣ`.
"""
function exp2(d::T) where {T<:Union{Float32,Float64}}
    q = round(d)
    qi = unsafe_trunc(Int, q)

    s = d - q

    u = exp2_kernel(s)
    u = T(dnormalize(dadd(T(1.0), dmul(u,s))))

    u = ldexp2k(u, qi)

    d > max_exp2(T) && (u = T(Inf))
    d < min_exp2(T) && (u = T(0.0))
    return u
end


const max_exp10(::Type{Float64}) = 3.08254715559916743851e2 # log 2^1023*(2-2^-52)
const max_exp10(::Type{Float32}) = 38.531839419103626f0 # log 2^127 *(2-2^-23) 

const min_exp10(::Type{Float64}) = -3.23607245338779784854769e2 # log10 2^-1075
const min_exp10(::Type{Float32}) = -45.15449934959718f0         # log10 2^-150

@inline function exp10_kernel(x::Float64)
    c11 = 0.2411463498334267652e-3
    c10 = 0.1157488415217187375e-2
    c9  = 0.5013975546789733659e-2
    c8  = 0.1959762320720533080e-1
    c7  = 0.6808936399446784138e-1
    c6  = 0.2069958494722676234e0
    c5  = 0.5393829292058536229e0
    c4  = 0.1171255148908541655e1
    c3  = 0.2034678592293432953e1
    c2  = 0.2650949055239205876e1
    c1  = 0.2302585092994045901e1
    return @horner x c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11
end

@inline function exp10_kernel(x::Float32)
    c6 = 0.2064004987f0
    c5 = 0.5417877436f0
    c4 = 0.1171286821f1
    c3 = 0.2034656048f1
    c2 = 0.2650948763f1
    c1 = 0.2302585125f1
    return @horner x c1 c2 c3 c4 c5 c6
end

"""
    exp10(x)

Compute the base-`10` exponential of `x`, that is `10ˣ`.
"""
function exp10(d::T) where {T<:Union{Float32,Float64}}
    q = round(T(MLOG10_2) * d)
    qi = unsafe_trunc(Int, q)

    s = muladd(q, -L10U(T), d)
    s = muladd(q, -L10L(T), s)

    u = exp10_kernel(s)
    u = T(dnormalize(dadd(T(1.0), dmul(u,s))))

    u = ldexp2k(u, qi)

    d > max_exp10(T) && (u = T(Inf))
    d < min_exp10(T) && (u = T(0.0))

    return u
end


const max_expm1(::Type{Float64}) = 7.09782712893383996732e2 # log 2^1023*(2-2^-52)
const max_expm1(::Type{Float32}) = 88.72283905206835f0 # log 2^127 *(2-2^-23)

const min_expm1(::Type{Float64}) = -37.42994775023704434602223
const min_expm1(::Type{Float32}) = -17.3286790847778338076068394f0

"""
    expm1(x)

Compute `eˣ- 1` accurately for small values of `x`.
"""
function expm1(x::T) where {T<:Union{Float32,Float64}}
    u = T(dadd2(expk2(Double(x)), -T(1.0)))
    x > max_expm1(T) && (u = T(Inf))
    x < min_expm1(T) && (u = -T(1.0))
    isnegzero(x) && (u = T(-0.0))
    return u
end


const max_exp(::Type{Float64}) = 709.78271114955742909217217426  # log 2^1023*(2-2^-52)
const max_exp(::Type{Float32}) = 88.72283905206835f0             # log 2^127 *(2-2^-23)

const min_exp(::Type{Float64}) = -7.451332191019412076235e2 # log 2^-1075
const min_exp(::Type{Float32}) = -103.97208f0               # ≈ log 2^-150

@inline function exp_kernel(x::Float64)
    c11 = 2.08860621107283687536341e-09
    c10 = 2.51112930892876518610661e-08
    c9  = 2.75573911234900471893338e-07
    c8  = 2.75572362911928827629423e-06
    c7  = 2.4801587159235472998791e-05
    c6  = 0.000198412698960509205564975
    c5  = 0.00138888888889774492207962
    c4  = 0.00833333333331652721664984
    c3  = 0.0416666666666665047591422
    c2  = 0.166666666666666851703837
    c1  = 0.50
    return @horner x c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11
end

@inline function exp_kernel(x::Float32)
    c6 = 0.000198527617612853646278381f0
    c5 = 0.00139304355252534151077271f0
    c4 = 0.00833336077630519866943359f0
    c3 = 0.0416664853692054748535156f0
    c2 = 0.166666671633720397949219f0
    c1 = 0.5f0
    return @horner x c1 c2 c3 c4 c5 c6
end

"""
    exp(x)

Compute the base-`e` exponential of `x`, that is `eˣ`.
"""
function exp(d::T) where {T<:Union{Float32,Float64}}
    q = round(T(MLN2E) * d)
    qi = unsafe_trunc(Int, q)

    s = muladd(q, -L2U(T), d)
    s = muladd(q, -L2L(T), s)

    u = exp_kernel(s)

    u = s * s * u + s + 1
    u = ldexp2k(u, qi)

    d > max_exp(T) && (u = T(Inf))
    d < min_exp(T) && (u = T(0))

    return u
end
