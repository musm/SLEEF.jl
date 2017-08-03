## utility functions mainly used by the private math functions in priv.jl

function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)), nextfloat(one(T)), -nextfloat(one(T), 2)) != zero(T))
end
const FMA_FAST = is_fma_fast(Float64) && is_fma_fast(Float32)

@inline isnegzero(x::T) where {T<:AbstractFloat} = x === T(-0.0)

@inline ispinf(x::T) where {T<:AbstractFloat} = x == T(Inf)
@inline isninf(x::T) where {T<:AbstractFloat} = x == T(-Inf)

# _sign emits better native code than sign but does not properly handle the Inf/NaN cases
@inline _sign(d::T) where {T<:AbstractFloat} = flipsign(one(T), d)

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline pow2i(::Type{T}, q::Int) where {T<:Union{Float32,Float64}} = integer2float(T, q + exponent_bias(T))

# sqrt without the domain checks which we don't need since we handle the checks ourselves
if VERSION < v"0.7-"
    _sqrt(x::T) where {T<:Union{Float32,Float64}} = Base.sqrt_llvm_fast(x)
else
    _sqrt(x::T) where {T<:Union{Float32,Float64}} = Base.sqrt_llvm(x)
end
