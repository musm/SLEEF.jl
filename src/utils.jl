## utility functions mainly used by the private math functions in priv.jl

function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)), nextfloat(one(T)), -nextfloat(one(T), 2)) != zero(T))
end
const FMA_FAST = is_fma_fast(Float64) && is_fma_fast(Float32)

@inline isnegzero(x::T) where {T<:Union{Float32,Float64}} = x === T(-0.0)
# Disabling the check for performance when vecterized.
# A PR succesfully vectorizing the check is welcome.
@inline isnegzero(x::Vec{N}) where {N} = Vec{N,Bool}(false)

@inline ispinf(x::T) where {T<:FloatType} = x == T(Inf)
@inline isninf(x::T) where {T<:FloatType} = x == T(-Inf)

# _sign emits better native code than sign but does not properly handle the Inf/NaN cases
@inline _sign(d::T) where {T<:FloatType} = flipsign(one(T), d)

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline function integer2float(::Type{<:Union{Vec{N,Float64},Float64}}, m::Vec{N,Int}) where {N}
    reinterpret(Vec{N,Float64}, Vec{N,Int64}(ntuple(Val(N)) do n
        Core.VecElement(m[n] % Int64)
    end) << significand_bits(Float64))
end
@inline function integer2float(::Type{<:Union{Vec{N,Float32},Float32}}, m::Vec{N,<:Integer}) where {N}
    reinterpret(Vec{N,Float32}, Vec{N,Int32}(ntuple(Val(N)) do n
        Core.VecElement(m[n] % Int32)
    end) << Int32(significand_bits(Float32)))
end

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline function float2integer(d::Vec{N,Float64}) where {N}
    Vec{N,Int64}(ntuple(Val(N)) do n
        Core.VecElement((reinterpret(Int64, d[n]) >> significand_bits(Float64)) % Int)
    end)
end
@inline function float2integer(d::Vec{N,Float32}) where {N}
    Vec{N,Int32}(ntuple(Val(N)) do n
        Core.VecElement((reinterpret(Int32, d[n]) >> Int32(significand_bits(Float32))) % Int32)
    end)
end

@inline pow2i(::Type{T}, q::I) where {T<:Union{Float32,Float64},I<:IntegerType} = integer2float(T, q + exponent_bias(T))

# sqrt without the domain checks which we don't need since we handle the checks ourselves
@inline _sqrt(x::T) where {T<:Union{Float32,Float64}} = Base.sqrt_llvm(x)
@inline _sqrt(x::Vec) = sqrt(x)
