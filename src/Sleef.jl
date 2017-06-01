__precompile__()

module SLEEF

# export sin, cos, tan, asin, acos, atan, atan2, sincos, sinh, cosh, tanh,
#     asinh, acosh, atanh, log, log2, log10, log1p, ilogb, exp, exp2, exp10, expm1, ldexp, cbrt, pow

# fast variants (within 3 ulp)
# export sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast

using Base.Math: @horner, exponent_bias, exponent_mask, significand_bits, IEEEFloat, exponent_raw_max

## constants

const MLN2  = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01 # log(2)
const MLN2E = 1.442695040888963407359924681001892137426645954152985934135449406931109219181187     # log2(e)

const MPI  = 3.141592653589793238462643383279502884197169399375105820974944592307816406286198     # pi
const MPI2 = 1.570796326794896619231321691639751442098584699687552910487472296153908203143099     # pi/2
const MPI4 = 7.853981633974483096156608458198757210492923498437764552437361480769541015715495e-01 # pi/4
const M1PI = 3.183098861837906715377675267450287240689192914809128974953346881177935952684543e-01 # 1/pi
const M2PI = 6.366197723675813430755350534900574481378385829618257949906693762355871905369086e-01 # 2/pi
const M4PI = 1.273239544735162686151070106980114896275677165923651589981338752471174381073817     # 4/pi

const MSQRT2 = 1.414213562373095048801688724209698078569671875376948073176679737990732478462102 # sqrt(2)
const M1SQRT2 = 7.071067811865475244008443621048490392848359376884740365883398689953662392310596e-01 # 1/sqrt(2)

const M2P13 = 1.259921049894873164767210607278228350570251464701507980081975112155299676513956 # 2^1/3
const M2P23 = 1.587401051968199474751705639272308260391493327899853009808285761825216505624206 # 2^2/3

MDLN10E(::Type{Float64}) = Double(0.4342944819032518, 1.098319650216765e-17) # log10(e)
MDLN10E(::Type{Float32}) = Double(0.4342945f0, -1.010305f-8)

MDLN2E(::Type{Float64}) = Double(1.4426950408889634, 2.0355273740931033e-17) # log2(e)
MDLN2E(::Type{Float32}) = Double(1.442695f0, 1.925963f-8)

MDLN10(::Type{Float64}) = Double(2.302585092994046, -2.1707562233822494e-16) # log(10)
MDLN10(::Type{Float32}) = Double(2.3025851f0, -3.1975436f-8)

MDLN2(::Type{Float64}) = Double(0.6931471805599453, 2.3190468138462996e-17)  # log(2)
MDLN2(::Type{Float32}) = Double(0.6931472f0, -1.9046542f-9)

MDPI(::Type{Float64})  = Double(3.141592653589793, 1.2246467991473532e-16) # pi
MDPI(::Type{Float32})  = Double(3.1415927f0, -8.742278f-8)
MDPI2(::Type{Float64}) = Double(1.5707963267948966, 6.123233995736766e-17) # pi/2
MDPI2(::Type{Float32}) = Double(1.5707964f0, -4.371139f-8)

MD2P13(::Type{Float64}) = Double(1.2599210498948732, -2.589933375300507e-17) # 2^1/3
MD2P13(::Type{Float32}) = Double(1.2599211f0, -2.4018702f-8)

MD2P23(::Type{Float64}) = Double(1.5874010519681996, -1.0869008194197823e-16) # 2^2/3
MD2P23(::Type{Float32}) = Double(1.587401f0, 1.9520385f-8)

# Split pi into four parts (each is 26 bits)
PIA(::Type{Float64}) = 3.1415926218032836914
PIB(::Type{Float64}) = 3.1786509424591713469e-08
PIC(::Type{Float64}) = 1.2246467864107188502e-16
PID(::Type{Float64}) = 1.2736634327021899816e-24

PIA(::Type{Float32}) = 3.140625f0
PIB(::Type{Float32}) = 0.0009670257568359375f0
PIC(::Type{Float32}) = 6.2771141529083251953f-7
PID(::Type{Float32}) = 1.2154201256553420762f-10

# split 2/pi into upper and lower parts
M2PIU(::Type{Float64}) = 0.63661977236758138243
M2PIL(::Type{Float64}) = -3.9357353350364971764e-17

# Split log(2) into upper and lower parts
LN2U(::Type{Float64}) = 0.69314718055966295651160180568695068359375
LN2L(::Type{Float64}) = 0.28235290563031577122588448175013436025525412068e-12

LN2U(::Type{Float32}) = 0.693145751953125f0
LN2L(::Type{Float32}) = 1.428606765330187045f-06

TRIG_MAX(::Type{Float64}) = 1e14
SQRT_MAX(::Type{Float64}) = 1.3407807929942596355e154

TRIG_MAX(::Type{Float32}) = 1f7
SQRT_MAX(::Type{Float32}) = 18446743523953729536f0

include("utils.jl")  # utility functions
include("double.jl") # Dekker style double double functions
include("priv.jl")   # private math functions
include("exp.jl")    # exponential functions
include("log.jl")    # logarithmic functions
include("trig.jl")   # trigonometric and inverse trigonometric functions
include("hyp.jl")    # hyperbolic and inverse hyperbolic functions
include("misc.jl")   # miscallenous math functions including pow and cbrt

# fallback definitions

for func in (:sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh,
             :asinh, :acosh, :atanh, :log, :log2, :log10, :log1p, :exp, :exp2, :exp10, :expm1, :cbrt,
             :sin_fast, :cos_fast, :tan_fast, :sincos_fast, :asin_fast, :acos_fast, :atan_fast, :atan2_fast, :log_fast, :cbrt_fast)
    @eval begin
        $func(a::Float16) = Float16($func(Float32(a)))
        $func(x::Real) = $func(float(x))
    end
end

for func in (:atan2, :hypot)
    @eval begin
        $func(y::Real, x::Real) = $func(promote(float(y), float(x))...)
        $func(a::Float16,b::Float16) = Float16($func(Float32(a),Float32(b)))
    end
end

end
