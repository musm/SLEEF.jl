@testset "exceptional $T" for T in (Float32, Float64)

@testset "exceptional $xatan" for xatan in (SLEEF.atan_fast, SLEEF.atan)

    @test xatan(T(0.0),  T(-0.0)) === T(pi)
    @test xatan(T(-0.0), T(-0.0)) === -T(pi)
    @test ispzero(xatan(T(0.0), T(0.0)))
    @test isnzero(xatan(T(-0.0), T(0.0)))
    @test xatan( T(Inf), -T(Inf)) === T(3*pi/4)
    @test xatan(-T(Inf), -T(Inf)) === T(-3*pi/4)
    @test xatan( T(Inf),  T(Inf))  === T(pi/4)
    @test xatan(-T(Inf),  T(Inf))  === T(-pi/4)


    y = T(0.0)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan(y,x) === T(pi)
    end


    y = T(-0.0)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan(y,x) === T(-pi)
    end


    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    xa = T[T(0.0), T(-0.0)]
    for x in xa, y in ya
        @test xatan(y,x) === T(-pi/2)
    end


    ya = T[100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5]
    xa = T[T(0.0), T(-0.0)]
    for x in xa, y in ya
        @test xatan(y,x) === T(pi/2)
    end


    y = T(Inf)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for x in xa
        @test xatan(y,x) === T(pi/2)
    end


    y = T(-Inf)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for x in xa
        @test xatan(y,x) === T(-pi/2)
    end


    ya = T[0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    x = T(Inf)
    for y in ya
        @test ispzero(xatan(y,x))
    end


    ya = T[-0.5, -1.5, -2.0, -2.5, -3.0, -100000, -100000.5]
    x = T(Inf)
    for y in ya
        @test isnzero(xatan(y,x))
    end


    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
    x = T(NaN)
    for y in ya
        @test isnan(xatan(y,x))
    end


    y = T(NaN)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
    for x in xa
        @test isnan(xatan(y,x))
    end

end # denormal/nonumber atan


@testset "exceptional xpow" begin

    @test SLEEF.pow(T(1), T(NaN))   === T(1)
    @test SLEEF.pow( T(NaN), T(0))  === T(1)
    @test SLEEF.pow(-T(1), T(Inf))  === T(1)
    @test SLEEF.pow(-T(1), T(-Inf)) === T(1)


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    ya = T[-100000.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 100000.5]
    for x in xa, y in ya
        @test isnan(SLEEF.pow(x,y))
    end


    x = T(NaN)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test isnan(SLEEF.pow(x,y))
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(NaN)
    for x in xa
        @test isnan(SLEEF.pow(x,y))
    end


    x = T(0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test ispzero(SLEEF.pow(x,y))
    end


    x = T(-0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test isnzero(SLEEF.pow(x,y))
    end


    xa = T[0.0, -0.0]
    ya = T[0.5, 1.5, 2.0, 2.5, 4.0, 100000, 100000.5]
    for x in xa, y in ya
        @test ispzero(SLEEF.pow(x,y))
    end


    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(-Inf)
    for x in xa
        @test SLEEF.pow(x,y) === T(Inf)
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(-Inf)
    for x in xa
        @test ispzero(SLEEF.pow(x,y))
    end


    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(Inf)
    for x in xa
        @test ispzero(SLEEF.pow(x,y))
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(Inf)
    for x in xa
        @test SLEEF.pow(x,y) === T(Inf)
    end


    x = T(-Inf)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test isnzero(SLEEF.pow(x,y))
    end


    x = T(-Inf)
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for y in ya
        @test ispzero(SLEEF.pow(x,y))
    end


    x = T(-Inf)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test SLEEF.pow(x,y) === T(-Inf)
    end


    x = T(-Inf)
    ya = T[0.5, 1.5, 2, 2.5, 3.5, 4, 100000, 100000.5]
    for y in ya
        @test SLEEF.pow(x,y) === T(Inf)
    end


    x = T(Inf)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for y in ya
        @test ispzero(SLEEF.pow(x,y))
    end


    x = T(Inf)
    ya = T[0.5, 1, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test SLEEF.pow(x,y) === T(Inf)
    end


    x = T(0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test SLEEF.pow(x,y) ===  T(Inf)
    end


    x = T(-0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test SLEEF.pow(x,y) ===  T(-Inf)
    end


    xa = T[0.0, -0.0]
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for x in xa, y in ya
        @test SLEEF.pow(x,y) === T(Inf)
    end


    xa = T[1000, -1000]
    ya = T[1000, 1000.5, 1001]
    for x in xa, y in ya
        @test cmpdenorm(SLEEF.pow(x,y), Base.:^(BigFloat(x),BigFloat(y)))
    end

end # denormal/nonumber pow


fun_table = Dict(SLEEF.sin_fast => Base.sin, SLEEF.sin => Base.sin)
@testset "exceptional $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


fun_table = Dict(SLEEF.cos_fast => Base.cos, SLEEF.cos => Base.cos)
@testset "exceptional $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


@testset "exceptional sin in $xsincos"for xsincos in (SLEEF.sincos_fast, SLEEF.sincos)
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        q = xsincos(x)[1]
        @test cmpdenorm(q, Base.sin(BigFloat(x)))
    end
end


@testset "exceptional cos in $xsincos"for xsincos in (SLEEF.sincos_fast, SLEEF.sincos)
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        q = xsincos(x)[2]
        @test cmpdenorm(q, Base.cos(BigFloat(x)))
    end
end


@testset "exceptional $xtan" for xtan in (SLEEF.tan_fast, SLEEF.tan)
    xa = T[NaN, Inf, -Inf, -0.0, 0.0, pi/2, -pi/2]
    for x in xa
        @test cmpdenorm(xtan(x), Base.tan(BigFloat(x)))
    end
end


fun_table = Dict(SLEEF.asin => Base.asin, SLEEF.asin_fast => Base.asin, SLEEF.acos => Base.acos, SLEEF.acos_fast => Base.acos)
@testset "exceptional $xatrig" for (xatrig, atrig) in fun_table
    xa = T[NaN, Inf, -Inf, 2, -2, 1, -1, -0.0, 0.0]
    for x in xa
        @test cmpdenorm(xatrig(x), atrig(BigFloat(x)))
    end
end


@testset "exceptional $xatan" for xatan in (SLEEF.atan, SLEEF.atan_fast)
    xa = T[NaN, Inf, -Inf, -0.0, 0.0]
    for x in xa
        @test cmpdenorm(xatan(x), Base.atan(BigFloat(x)))
    end
end


@testset "exceptional exp" begin
    xa = T[NaN, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.exp(x), Base.exp(BigFloat(x)))
    end
end


@testset "exceptional sinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.sinh(x), Base.sinh(BigFloat(x)))
    end
end


@testset "exceptional cosh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.cosh(x),  Base.cosh(BigFloat(x)))
    end
end


@testset "exceptional tanh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.tanh(x),  Base.tanh(BigFloat(x)))
    end
end


@testset "exceptional asinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.asinh(x), Base.asinh(BigFloat(x)))
    end
end


@testset "exceptional acosh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.acosh(x), Base.acosh(BigFloat(x)))
    end
end


@testset "exceptional atanh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, -1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(SLEEF.atanh(x), Base.atanh(BigFloat(x)))
    end
end


@testset "exceptional $xcbrt" for xcbrt = (SLEEF.cbrt, SLEEF.cbrt_fast)
    xa = T[NaN, Inf, -Inf, 0.0, -0.0]
    for x in xa
        @test cmpdenorm(SLEEF.cbrt(x), Base.cbrt(BigFloat(x)))
    end
end


@testset "exceptional exp2" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(SLEEF.exp2(x), Base.exp2(BigFloat(x)))
    end
end


@testset "exceptional exp10" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(SLEEF.exp10(x), Base.exp10(BigFloat(x)))
    end
end


@testset "exceptional expm1" begin
    xa = T[NaN, Inf, -Inf, 0.0, -0.0]
    for x in xa
        @test cmpdenorm(SLEEF.expm1(x), Base.expm1(BigFloat(x)))
    end
end



@testset "exceptional $xlog" for xlog in (SLEEF.log, SLEEF.log_fast)
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(xlog(x), Base.log(BigFloat(x)))
    end
end


@testset "exceptional log10" begin
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(SLEEF.log10(x), Base.log10(BigFloat(x)))
    end
end


@testset "exceptional log2" begin
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(SLEEF.log2(x), Base.log2(BigFloat(x)))
    end
end


@testset "exceptional log1p" begin
    xa = T[NaN, Inf, -Inf, 0.0, -0.0, -1.0, -2.0]
    for x in xa
        @test cmpdenorm(SLEEF.log1p(x), Base.log1p(BigFloat(x)))
    end
end


@testset "exceptional ldexp" begin
    for i = -10000:10000
        a = SLEEF.ldexp(T(1.0), i)
        b = Base.ldexp(BigFloat(1.0), i)
        @test (isfinite(b) && a == b || cmpdenorm(a,b))
    end
end


@testset "exceptional ilogb" begin
    @test SLEEF.ilogb(+T(Inf)) == SLEEF.INT_MAX
    @test SLEEF.ilogb(-T(Inf)) == SLEEF.INT_MAX
    @test SLEEF.ilogb(+T(0.0)) == SLEEF.FP_ILOGB0
    @test SLEEF.ilogb(-T(0.0)) == SLEEF.FP_ILOGB0
    @test SLEEF.ilogb( T(NaN)) == SLEEF.FP_ILOGBNAN
end



end #exceptional
