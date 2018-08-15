using SLEEF
using BenchmarkTools
using JLD, DataStructures
using Printf

const RETUNE  = false
const VERBOSE = true
const DETAILS = false
const test_types = (Float64, Float32) # Which types do you want to bench?
const bench = ("Base", "SLEEF")
const suite = BenchmarkGroup()
for n in bench
    suite[n] = BenchmarkGroup([n])
end

bench_reduce(f::Function, X) = mapreduce(x -> reinterpret(Unsigned,x), |, f(x) for x in X)

using Base.Math.IEEEFloat

MRANGE(::Type{Float64}) = 10000000
MRANGE(::Type{Float32}) = 10000
IntF(::Type{Float64}) = Int64
IntF(::Type{Float32}) = Int32
x_trig(::Type{T}) where {T<:IEEEFloat} = begin
    x_trig = T[]
    for i = 1:10000
        s = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) - IntF(T)(20))
        e = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) + IntF(T)(20))
        d = s
        while d <= e 
            append!(x_trig, d)
            d = reinterpret(T, reinterpret(IntF(T), d) + IntF(T)(1))
        end
    end
    x_trig = append!(x_trig, -10:0.0002:10)
    x_trig = append!(x_trig, -MRANGE(T):200.1:MRANGE(T))
end
x_exp(::Type{T}) where {T<:IEEEFloat}        = map(T, vcat(-10:0.0002:10, -1000:0.1:1000))
x_exp2(::Type{T}) where {T<:IEEEFloat}       = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
x_exp10(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
x_expm1(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(-10:0.0002:10, -1000:0.021:1000, -1000:0.023:1000, 10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300), 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300)))
x_log(::Type{T}) where {T<:IEEEFloat}        = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
x_log10(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000))
x_log1p(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000, 10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300)))
x_atrig(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(-1:0.00002:1))
x_atan(::Type{T}) where {T<:IEEEFloat}       = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
x_cbrt(::Type{T}) where {T<:IEEEFloat}       = map(T, vcat(-10000:0.2:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
x_trigh(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(-10:0.0002:10, -1000:0.02:1000))
x_asinhatanh(::Type{T}) where {T<:IEEEFloat} = map(T, vcat(-10:0.0002:10, -1000:0.02:1000))
x_acosh(::Type{T}) where {T<:IEEEFloat}      = map(T, vcat(1:0.0002:10, 1:0.02:1000))
x_pow(::Type{T}) where {T<:IEEEFloat} = begin
    xx1 = map(Tuple{T,T}, [(x,y) for x = -100:0.20:100, y = 0.1:0.20:100])[:]
    xx2 = map(Tuple{T,T}, [(x,y) for x = -100:0.21:100, y = 0.1:0.22:100])[:]
    xx3 = map(Tuple{T,T}, [(x,y) for x = 2.1, y = -1000:0.1:1000])
    xx = vcat(xx1, xx2, xx2)
end

import Base.atanh
for f in (:atanh,)
    @eval begin
        ($f)(x::Float64) = ccall($(string(f)), Float64, (Float64,), x)
        ($f)(x::Float32) = ccall($(string(f,"f")), Float32, (Float32,), x)
    end
end

const micros = OrderedDict(
    "sin"   => x_trig,
    "cos"   => x_trig,
    "tan"   => x_trig,
    "asin"  => x_atrig,
    "acos"  => x_atrig,
    "atan"  => x_atan,
    "exp"   => x_exp,
    "exp2"  => x_exp2,
    "exp10" => x_exp10,
    "expm1" => x_expm1,
    "log"   => x_log,
    "log2"  => x_log10,
    "log10" => x_log10,
    "log1p" => x_log1p,
    "sinh"  => x_trigh,
    "cosh"  => x_trigh,
    "tanh"  => x_trigh,
    "asinh" => x_asinhatanh,
    "acosh" => x_acosh,
    "atanh" => x_asinhatanh,
    "cbrt"  => x_cbrt
    )

for n in bench
    for (f,x) in micros
        suite[n][f] = BenchmarkGroup([f])
        for T in test_types
            fex = Expr(:., Symbol(n), QuoteNode(Symbol(f)))
            suite[n][f][string(T)] = @benchmarkable bench_reduce($fex, $(x(T)))
        end
    end
end


tune_params = joinpath(@__DIR__, "params.jld")
if !isfile(tune_params) || RETUNE
    tune!(suite; verbose=VERBOSE, seconds = 2)
    save(tune_params, "suite", params(suite))
    println("Saving tuned parameters.")
else
    println("Loading pretuned parameters.")
    loadparams!(suite, load(tune_params, "suite"), :evals, :samples)
end

println("Running micro benchmarks...")
results = run(suite; verbose=VERBOSE, seconds = 2)

printstyled("Benchmarks: median ratio SLEEF/Base\n", color = :blue)
for f in keys(micros)
    printstyled(string(f) color = :magenta)
    for T in test_types
        println()
        print("time: ", )
        tratio = ratio(median(results["SLEEF"][f][string(T)]), median(results["Base"][f][string(T)])).time
        tcolor = tratio > 3 ? :red : tratio < 1.5 ? :green : :blue
        printstyled(@sprintf("%.2f",tratio), " ", string(T), color = tcolor)
        if DETAILS
            printstyled("details SLEEF/Base\n", color=:blue)
            println(results["SLEEF"][f][string(T)])
            println(results["Base"][f][string(T)])
            println()
        end
    end
    println("\n")
end
