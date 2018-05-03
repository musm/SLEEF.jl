import Base: -, <, copysign, flipsign, convert

struct Double{T<:IEEEFloat} <: Number
    hi::T
    lo::T
end
Double(x::T) where {T<:IEEEFloat} = Double(x, zero(T))

(::Type{T})(x::Double{T}) where {T<:IEEEFloat} = x.hi + x.lo


@inline trunclo(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xffff_ffff_f800_0000) # clear lower 27 bits (leave upper 26 bits)
@inline trunclo(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xffff_f000) # clear lowest 12 bits (leave upper 12 bits)

@inline function splitprec(x::IEEEFloat)
    hx = trunclo(x)
    hx, x - hx
end

@inline function dnormalize(x::Double{T}) where {T}
    r = x.hi + x.lo
    Double(r, (x.hi - r) + x.lo)
end

@inline flipsign(x::Double{T}, y::T) where {T<:IEEEFloat} = Double(flipsign(x.hi, y), flipsign(x.lo, y))

@inline scale(x::Double{T}, s::T) where {T<:IEEEFloat} = Double(s * x.hi, s * x.lo)


@inline (-)(x::Double{T}) where {T<:IEEEFloat} = Double(-x.hi, -x.lo)

@inline function (<)(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
    x.hi < y.hi
end

@inline function (<)(x::Double{T}, y::Number) where {T<:IEEEFloat}
    x.hi < y
end

@inline function (<)(x::Number, y::Double{T}) where {T<:IEEEFloat}
    x < y.hi
end

# quick-two-sum x+y
@inline function dadd(x::T, y::T) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x + y
    Double(s, (x - s) + y)
end

@inline function dadd(x::T, y::Double{T}) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x + y.hi
    Double(s, (x - s) + y.hi + y.lo)
end

@inline function dadd(x::Double{T}, y::T) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x.hi + y
    Double(s, (x.hi - s) + y + x.lo)
end

@inline function dadd(x::Double{T}, y::Double{T}) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x.hi + y.hi
    Double(s, (x.hi - s) + y.hi + y.lo + x.lo)
end

@inline function dsub(x::Double{T}, y::Double{T}) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x.hi - y.hi
    Double(s, (x.hi - s) - y.hi - y.lo + x.lo)
end

@inline function dsub(x::Double{T}, y::T) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x.hi - y
    Double(s, (x.hi - s) - y + x.lo)
end

@inline function dsub(x::T, y::Double{T}) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x - y.hi
    Double(s, (x - s) - y.hi - y.lo)
end

@inline function dsub(x::T, y::T) where {T<:IEEEFloat} #WARNING |x| >= |y|
    s = x - y
    Double(s, (x - s) - y)
end


# two-sum x+y  NO BRANCH
@inline function dadd2(x::T, y::T) where {T<:IEEEFloat}
    s = x + y
    v = s - x
    Double(s, (x - (s - v)) + (y - v))
end

@inline function dadd2(x::T, y::Double{T}) where {T<:IEEEFloat}
    s = x + y.hi
    v = s - x
    Double(s, (x - (s - v)) + (y.hi - v) + y.lo)
end

@inline dadd2(x::Double{T}, y::T) where {T<:IEEEFloat} = dadd2(y, x)

@inline function dadd2(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
    s = x.hi + y.hi
    v = s - x.hi
    Double(s, (x.hi - (s - v)) + (y.hi - v) + x.lo + y.lo)
end

@inline function dsub2(x::T, y::T) where {T<:IEEEFloat}
    s = x - y
    v = s - x
    Double(s, (x - (s - v)) + (-y - v))
end

@inline function dsub2(x::T, y::Double{T}) where {T<:IEEEFloat}
    s = x - y.hi
    v = s - x
    Double(s, (x - (s - v)) + (-y.hi - v) - y.lo)
end

@inline function dsub2(x::Double{T}, y::T) where {T<:IEEEFloat}
    s = x.hi - y
    v = s - x.hi
    Double(s, (x.hi - (s - v)) + (-y - v) + x.lo)
end

@inline function dsub2(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
    s = x.hi - y.hi
    v = s - x.hi
    Double(s, (x.hi - (s - v)) + (-y.hi - v) + x.lo - y.lo)
end



if FMA_FAST

    # two-prod-fma
    @inline function dmul(x::T, y::T) where {T<:IEEEFloat}
        z = x * y
        Double(z, fma(x, y, -z))
    end

    @inline function dmul(x::Double{T}, y::T) where {T<:IEEEFloat}
        z = x.hi * y
        Double(z, fma(x.hi, y, -z) + x.lo * y)
    end

    @inline dmul(x::T, y::Double{T}) where {T<:IEEEFloat} = dmul(y, x)

    @inline function dmul(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
        z = x.hi * y.hi
        Double(z, fma(x.hi, y.hi, -z) + x.hi * y.lo + x.lo * y.hi)
    end

    # x^2
    @inline function dsqu(x::T) where {T<:IEEEFloat}
        z = x * x
        Double(z, fma(x, x, -z))
    end

    @inline function dsqu(x::Double{T}) where {T<:IEEEFloat}
        z = x.hi * x.hi
        Double(z, fma(x.hi, x.hi, -z) + x.hi * (x.lo + x.lo))
    end

    # sqrt(x)
    @inline function dsqrt(x::Double{T}) where {T<:IEEEFloat}
        zhi = _sqrt(x.hi)
        Double(zhi, (x.lo + fma(-zhi, zhi, x.hi)) / (zhi + zhi))
    end

    # x/y
    @inline function ddiv(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
        invy = 1 / y.hi
        zhi = x.hi * invy
        Double(zhi, (fma(-zhi, y.hi, x.hi) + fma(-zhi, y.lo, x.lo)) * invy)
    end

    @inline function ddiv(x::T, y::T) where {T<:IEEEFloat}
        ry = 1 / y
        r = x * ry
        Double(r, fma(-r, y, x) * ry)
    end

    # 1/x
    @inline function drec(x::T) where {T<:IEEEFloat}
        zhi = 1 / x
        Double(zhi, fma(-zhi, x, one(T)) * zhi)
    end

    @inline function drec(x::Double{T}) where {T<:IEEEFloat}
        zhi = 1 / x.hi
        Double(zhi, (fma(-zhi, x.hi, one(T)) + -zhi * x.lo) * zhi)
    end

else

    #two-prod x*y
    @inline function dmul(x::T, y::T) where {T<:IEEEFloat}
        hx, lx = splitprec(x)
        hy, ly = splitprec(y)
        z = x * y
        Double(z, ((hx * hy - z) + lx * hy + hx * ly) + lx * ly)
    end

    @inline function dmul(x::Double{T}, y::T) where {T<:IEEEFloat}
        hx, lx = splitprec(x.hi)
        hy, ly = splitprec(y)
        z = x.hi * y
        Double(z, (hx * hy - z) + lx * hy + hx * ly + lx * ly + x.lo * y)
    end

    @inline dmul(x::T, y::Double{T}) where {T<:IEEEFloat} = dmul(y, x)

    @inline function dmul(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
        hx, lx = splitprec(x.hi)
        hy, ly = splitprec(y.hi)
        z = x.hi * y.hi
        Double(z, (((hx * hy - z) + lx * hy + hx * ly) + lx * ly) + x.hi * y.lo + x.lo * y.hi)
    end

    # x^2
    @inline function dsqu(x::T) where {T<:IEEEFloat}
        hx, lx = splitprec(x)
        z = x * x
        Double(z, (hx * hx - z) + lx * (hx + hx) + lx * lx)
    end

    @inline function dsqu(x::Double{T}) where {T<:IEEEFloat}
        hx, lx = splitprec(x.hi)
        z = x.hi * x.hi
        Double(z, (hx * hx - z) + lx * (hx + hx) + lx * lx + x.hi * (x.lo + x.lo))
    end

    # sqrt(x)
    @inline function dsqrt(x::Double{T}) where {T<:IEEEFloat}
        c = _sqrt(x.hi)
        u = dsqu(c)
        Double(c, (x.hi - u.hi - u.lo + x.lo) / (c + c))
    end

    # x/y
    @inline function ddiv(x::Double{T}, y::Double{T}) where {T<:IEEEFloat}
        invy = 1 / y.hi
        c = x.hi * invy
        u = dmul(c, y.hi)
        Double(c, ((((x.hi - u.hi) - u.lo) + x.lo) - c * y.lo) * invy)
    end

    @inline function ddiv(x::T, y::T) where {T<:IEEEFloat}
        ry = 1 / y
        r = x * ry
        hx, lx = splitprec(r)
        hy, ly = splitprec(y)
        Double(r, (((-hx * hy + r * y) - lx * hy - hx * ly) - lx * ly) * ry)
    end


    # 1/x
    @inline function drec(x::T) where {T<:IEEEFloat}
        c = 1 / x
        u = dmul(c, x)
        Double(c, (one(T) - u.hi - u.lo) * c)
    end

    @inline function drec(x::Double{T}) where {T<:IEEEFloat}
        c = 1 / x.hi
        u = dmul(c, x.hi)
        Double(c, (one(T) - u.hi - u.lo - c * x.lo) * c)
    end

end
