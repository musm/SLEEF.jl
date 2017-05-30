<div align="center"> <img
src="https://rawgit.com/musm/Sleef.jl/master/doc/src/assets/logo.svg"
alt="SLEEF Logo" width="420"></img> </div>


A pure Julia port of the [SLEEF math library](https://github.com/shibatch/sleef).

**History**
- Release [v0.2.0](https://github.com/musm/Sleef.jl/releases/tag/v0.2.0) based on SLEEF v2.90
- Release [v0.1.0](https://github.com/musm/Sleef.jl/releases/tag/v0.1.0) based on SLEEF v2.80

<br><br>
[![Travis Build Status](https://travis-ci.org/musm/Sleef.jl.svg?branch=master)](https://travis-ci.org/musm/Sleef.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/j7lpafn4uf1trlfi/branch/master?svg=true)](https://ci.appveyor.com/project/musm/sleef-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/musm/Sleef.jl/badge.svg?branch=master)](https://coveralls.io/github/musm/Sleef.jl?branch=master)
[![codecov.io](http://codecov.io/github/musm/Sleef.jl/coverage.svg?branch=master)](http://codecov.io/github/musm/Sleef.jl?branch=master)

# Usage

If you have not already done so, we recommend building a custom system image by running
```julia
julia> if is_windows(); Pkg.add("WinRPM"); WinRPM.install("gcc"); end
julia> include(joinpath(dirname(JULIA_HOME), "share", "julia", "build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then restart `julia`. This will ensure you are taking full advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set)  if your CPU supports it.


To use  `Sleef.jl`
```julia
julia> Pkg.clone("https://github.com/JuliaMath/Sleef.jl.git")
julia> using Sleef

julia> Sleef.exp(3.0)
20.085536923187668

julia> Sleef.exp(3f0)
20.085537f0
```

The available functions include (within 1 ulp)
```julia
sin, cos, tan, asin, acos, atan, atan2, sincos, sinh, cosh, tanh,
    asinh, acosh, atanh, log, log2, log10, log1p, ilogb, exp, exp2, exp10, expm1, ldexp, cbrt, pow
 ```

Faster variants (within 3 ulp)
 ```julia
sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast
```

## Notes

The trigonometric functions are tested to return values with specified
accuracy when the argument is within the following range:

- Double (Float64) precision trigonometric functions : `|arg| <= 10000000`
- Single (Float32) precision trigonometric functions : `|arg| <= 10000`

# Benchmarks

You can benchmark the performance of the Sleef.jl math library on your machine by running

```julia
include(joinpath(Pkg.dir("Sleef"), "benchmark", "benchmark.jl"))
```

