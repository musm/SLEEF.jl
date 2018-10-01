<div align="center"> <img
src="https://rawgit.com/musm/SLEEF.jl/master/doc/src/assets/logo.svg"
alt="SLEEF Logo" width="380"></img> </div>


A pure Julia port of the [SLEEF math library](https://github.com/shibatch/SLEEF)

**History**
- Release [v0.4.0](https://github.com/musm/SLEEF.jl/releases/tag/v0.4.0) based on SLEEF v2.110
- Release [v0.3.0](https://github.com/musm/SLEEF.jl/releases/tag/v0.3.0) based on SLEEF v2.100
- Release [v0.2.0](https://github.com/musm/SLEEF.jl/releases/tag/v0.2.0) based on SLEEF v2.90
- Release [v0.1.0](https://github.com/musm/SLEEF.jl/releases/tag/v0.1.0) based on SLEEF v2.80

<br><br>
[![Travis Build Status](https://travis-ci.org/musm/SLEEF.jl.svg?branch=master)](https://travis-ci.org/musm/SLEEF.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/j7lpafn4uf1trlfi/branch/master?svg=true)](https://ci.appveyor.com/project/musm/SLEEF-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/musm/SLEEF.jl/badge.svg?branch=master)](https://coveralls.io/github/musm/SLEEF.jl?branch=master)
[![codecov.io](http://codecov.io/github/musm/SLEEF.jl/coverage.svg?branch=master)](http://codecov.io/github/musm/SLEEF.jl?branch=master)

# Usage

To use  `SLEEF.jl`
```julia
pkg> add SLEEF
julia> using SLEEF

julia> SLEEF.exp(3.0)
20.085536923187668

julia> SLEEF.exp(3f0)
20.085537f0
```

The available functions include (within 1 ulp)
```julia
sin, cos, tan, asin, acos, atan, sincos, sinh, cosh, tanh,
    asinh, acosh, atanh, log, log2, log10, log1p, ilogb, exp, exp2, exp10, expm1, ldexp, cbrt, pow
 ```

Faster variants (within 3 ulp)
 ```julia
sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast
```

## Notes

The trigonometric functions are tested to return values with specified
accuracy when the argument is within the following range:

- Double (Float64) precision trigonometric functions : `[-1e+14, 1e+14]`
- Single (Float32) precision trigonometric functions : `[-39000, 39000]`
