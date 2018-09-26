module SLEEF

# export sin, cos, tan, asin, acos, atan, sincos, sinh, cosh, tanh,
#     asinh, acosh, atanh, log, log2, log10, log1p, ilogb, exp, exp2, exp10, expm1, ldexp, cbrt, pow

# fast variants (within 3 ulp)
# export sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast

using Base.Math: uinttype, @horner, exponent_bias, exponent_mask, significand_bits, IEEEFloat, exponent_raw_max

## constants

const MLN2  = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01 # log(2)
const MLN2E = 1.442695040888963407359924681001892137426645954152985934135449406931 # log2(e)

const M_PI  = 3.141592653589793238462643383279502884 # pi
const PI_2 = 1.570796326794896619231321691639751442098584699687552910487472296153908203143099     # pi/2
const PI_4 = 7.853981633974483096156608458198757210492923498437764552437361480769541015715495e-01 # pi/4
const M_1_PI = 0.318309886183790671537767526745028724 # 1/pi
const M_2_PI = 0.636619772367581343075535053490057448 # 2/pi
const M_4_PI = 1.273239544735162686151070106980114896275677165923651589981338752471174381073817     # 4/pi

const MSQRT2 = 1.414213562373095048801688724209698078569671875376948073176679737990732478462102 # sqrt(2)
const M1SQRT2 = 7.071067811865475244008443621048490392848359376884740365883398689953662392310596e-01 # 1/sqrt(2)

const M2P13 = 1.259921049894873164767210607278228350570251464701507980081975112155299676513956 # 2^1/3
const M2P23 = 1.587401051968199474751705639272308260391493327899853009808285761825216505624206 # 2^2/3

const MLOG10_2 = 3.3219280948873623478703194294893901758648313930

const MDLN10E(::Type{Float64}) = Double(0.4342944819032518, 1.098319650216765e-17) # log10(e)
const MDLN10E(::Type{Float32}) = Double(0.4342945f0, -1.010305f-8)

const MDLN2E(::Type{Float64}) = Double(1.4426950408889634, 2.0355273740931033e-17) # log2(e)
const MDLN2E(::Type{Float32}) = Double(1.442695f0, 1.925963f-8)

const MDLN2(::Type{Float64}) = Double(0.693147180559945286226764, 2.319046813846299558417771e-17)  # log(2)
const MDLN2(::Type{Float32}) = Double(0.69314718246459960938f0, -1.904654323148236017f-9)

const MDPI(::Type{Float64})  = Double(3.141592653589793, 1.2246467991473532e-16) # pi
const MDPI(::Type{Float32})  = Double(3.1415927f0, -8.742278f-8)
const MDPI2(::Type{Float64}) = Double(1.5707963267948966, 6.123233995736766e-17) # pi/2
const MDPI2(::Type{Float32}) = Double(1.5707964f0, -4.371139f-8)

const MD2P13(::Type{Float64}) = Double(1.2599210498948732, -2.589933375300507e-17) # 2^1/3
const MD2P13(::Type{Float32}) = Double(1.2599211f0, -2.4018702f-8)

const MD2P23(::Type{Float64}) = Double(1.5874010519681996, -1.0869008194197823e-16) # 2^2/3
const MD2P23(::Type{Float32}) = Double(1.587401f0, 1.9520385f-8)

# Split pi into four parts (each is 26 bits)
const PI_A(::Type{Float64}) = 3.1415926218032836914
const PI_B(::Type{Float64}) = 3.1786509424591713469e-08
const PI_C(::Type{Float64}) = 1.2246467864107188502e-16
const PI_D(::Type{Float64}) = 1.2736634327021899816e-24

const PI_A(::Type{Float32}) = 3.140625f0
const PI_B(::Type{Float32}) = 0.0009670257568359375f0
const PI_C(::Type{Float32}) = 6.2771141529083251953f-7
const PI_D(::Type{Float32}) = 1.2154201256553420762f-10

const PI_XD(::Type{Float32}) = 1.2141754268668591976f-10
const PI_XE(::Type{Float32}) = 1.2446743939339977025f-13

# split 2/pi into upper and lower parts
const M_2_PI_H = 0.63661977236758138243
const M_2_PI_L = -3.9357353350364971764e-17

# Split log(10) into upper and lower parts
const L10U(::Type{Float64}) = 0.30102999566383914498
const L10L(::Type{Float64}) = 1.4205023227266099418e-13

const L10U(::Type{Float32}) = 0.3010253906f0
const L10L(::Type{Float32}) = 4.605038981f-6

# Split log(2) into upper and lower parts
const L2U(::Type{Float64}) = 0.69314718055966295651160180568695068359375
const L2L(::Type{Float64}) = 0.28235290563031577122588448175013436025525412068e-12

const L2U(::Type{Float32}) = 0.693145751953125f0
const L2L(::Type{Float32}) = 1.428606765330187045f-06

const TRIG_MAX(::Type{Float64}) = 1e14
const TRIG_MAX(::Type{Float32}) = 1f7

const SQRT_MAX(::Type{Float64}) = 1.3407807929942596355e154
const SQRT_MAX(::Type{Float32}) = 18446743523953729536f0

include("utils.jl")  # utility functions
include("double.jl") # Dekker style double double functions
include("priv.jl")   # private math functions
include("exp.jl")    # exponential functions
include("log.jl")    # logarithmic functions
include("trig.jl")   # trigonometric and inverse trigonometric functions
include("hyp.jl")    # hyperbolic and inverse hyperbolic functions
include("misc.jl")   # miscallenous math functions including pow and cbrt

# fallback definitions

for func in (:sin, :cos, :tan, :sincos, :asin, :acos, :atan, :sinh, :cosh, :tanh,
             :asinh, :acosh, :atanh, :log, :log2, :log10, :log1p, :exp, :exp2, :exp10, :expm1, :cbrt,
             :sin_fast, :cos_fast, :tan_fast, :sincos_fast, :asin_fast, :acos_fast, :atan_fast, :atan2_fast, :log_fast, :cbrt_fast)
    @eval begin
        $func(a::Float16) = Float16.($func(Float32(a)))
        $func(x::Real) = $func(float(x))
    end
end

for func in (:atan, :hypot)
    @eval begin
        $func(y::Real, x::Real) = $func(promote(float(y), float(x))...)
        $func(a::Float16, b::Float16) = Float16($func(Float32(a), Float32(b)))
    end
end
ldexp(x::Float16, q::Int) = Float16(ldexpk(Float32(x), q))

end
