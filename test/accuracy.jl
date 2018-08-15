
MRANGE(::Type{Float64}) = 10000000
MRANGE(::Type{Float32}) = 10000
IntF(::Type{Float64}) = Int64
IntF(::Type{Float32}) = Int32


@testset "Accuracy (max error in ulp) for $T" for T in (Float32, Float64)
    println("Accuracy tests for $T")


    xx = map(T, vcat(-10:0.0002:10, -1000:0.1:1000))
    fun_table = Dict(SLEEF.exp => Base.exp)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -1000:0.02:1000))
    fun_table = Dict(SLEEF.asinh => Base.asinh, SLEEF.atanh => Base.atanh)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(1:0.0002:10, 1:0.02:1000))
    fun_table = Dict(SLEEF.acosh => Base.acosh)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = T[]
    for i = 1:10000
        s = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) - IntF(T)(20))
        e = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) + IntF(T)(20))
        d = s
        while d <= e 
            append!(xx, d)
            d = reinterpret(T, reinterpret(IntF(T), d) + IntF(T)(1))
        end
    end
    xx = append!(xx, -10:0.0002:10)
    xx = append!(xx, -MRANGE(T):200.1:MRANGE(T))

    fun_table = Dict(SLEEF.sin => Base.sin, SLEEF.cos => Base.cos, SLEEF.tan => Base.tan)
    tol = 1
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.sin_fast => Base.sin, SLEEF.cos_fast => Base.cos, SLEEF.tan_fast => Base.tan)
    tol = 4
    test_acc(T, fun_table, xx, tol)

    global sin_sincos_fast(x) = (SLEEF.sincos_fast(x))[1]
    global cos_sincos_fast(x) = (SLEEF.sincos_fast(x))[2]
    fun_table = Dict(sin_sincos_fast => Base.sin, cos_sincos_fast => Base.cos)
    tol = 4
    test_acc(T, fun_table, xx, tol) 

    global sin_sincos(x) = (SLEEF.sincos(x))[1]
    global cos_sincos(x) = (SLEEF.sincos(x))[2]
    fun_table = Dict(sin_sincos => Base.sin, cos_sincos => Base.cos)
    tol = 1
    test_acc(T, fun_table, xx, tol) 
    
    
    xx = map(T, vcat(-1:0.00002:1))
    fun_table = Dict(SLEEF.asin_fast => Base.asin, SLEEF.acos_fast => Base.acos)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.asin => asin, SLEEF.acos => Base.acos)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
    fun_table = Dict(SLEEF.atan_fast => Base.atan)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.atan => Base.atan)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx1 = map(Tuple{T,T}, [zip(-10:0.050:10, -10:0.050:10)...])
    xx2 = map(Tuple{T,T}, [zip(-10:0.051:10, -10:0.052:10)...])
    xx3 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.51:100)...])
    xx4 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.52:100)...])
    xx = vcat(xx1, xx2, xx3, xx4)

    fun_table = Dict(SLEEF.atan_fast => Base.atan)
    tol = 2.5
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.atan => Base.atan)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    fun_table = Dict(SLEEF.log_fast => Base.log)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.log => Base.log)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000))
    fun_table = Dict(SLEEF.log10 => Base.log10, SLEEF.log2 => Base.log2)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000, 10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300)))
    fun_table = Dict(SLEEF.log1p => Base.log1p)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx1 = map(Tuple{T,T}, [(x,y) for x = -100:0.20:100, y = 0.1:0.20:100])[:]
    xx2 = map(Tuple{T,T}, [(x,y) for x = -100:0.21:100, y = 0.1:0.22:100])[:]
    xx3 = map(Tuple{T,T}, [(x,y) for x = 2.1, y = -1000:0.1:1000])
    xx = vcat(xx1, xx2, xx2)
    fun_table = Dict(SLEEF.pow => Base.:^)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10000:0.2:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    fun_table = Dict(SLEEF.cbrt_fast => Base.cbrt)
    tol = 2
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(SLEEF.cbrt => Base.cbrt)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
    fun_table = Dict(SLEEF.exp2 => Base.exp2)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
    fun_table = Dict(SLEEF.exp10 => Base.exp10)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -1000:0.021:1000, -1000:0.023:1000, 
        10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300), 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300)))
    fun_table = Dict(SLEEF.expm1 => Base.expm1)
    tol = 2
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -1000:0.02:1000))   
    fun_table = Dict(SLEEF.sinh => Base.sinh, SLEEF.cosh => Base.cosh, SLEEF.tanh => Base.tanh)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    @testset "xilogb at arbitrary values" begin
        xd = Dict{T,Int}(T(1e-30) => -100, T(2.31e-11) => -36, T(-1.0) => 0, T(1.0) => 0, 
                T(2.31e11) => 37,  T(1e30) => 99)
        for (i,j) in xd
            @test SLEEF.ilogb(i)  === j
        end
    end
    
end
