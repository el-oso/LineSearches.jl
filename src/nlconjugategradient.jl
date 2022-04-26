
# Linesearch
# LineSearch(this, reltol, abstol, initialbracketingrange)



# 1. Calculate steepeste descent
# 2. Set initial search direction to steepest descent -> g0
# 3. Set d0 <- g0
# 4. 

const Φ = Base.MathConstants.φ # Golden Ratio
const ϵ_min = 1e-21 # Tolerance to avoid division by zero

# Need a better place for this
const growLimit = 100
const maxEvaluations = 500

# function bracket_minimum(f, x=0; s=1e-2, k=2.0)
#     a, ya = x, f(x)
#     b, yb = a+s, f(a+s)
#     if yb > ya
#         a, b = b, a
#         ya, yb = yb, ya
#         s = -s
#     end

#     while true
#         c, yc = b+s,  f(b+s)
#         if yc > yb
#             return a < c? (a,c) : (c, a)
#         end
#         a, ya, b, yb = b, yb, c, yc
#         s *= k
#     end
# end

bisection(df, a::Number, b::Number)= bisection(df, promote(a, b)...)
bisection(df, a, b, ϵ) = bisection(df, a, b, ϵ=ϵ)
function bisection(df, a::T, b::T; ϵ=eps(T)) where T <: Number
    if a > b; a,b = b, a; end # ensure a < b

    ya, yb = df(a), df(b)
    if ya == zero(T); b = a; end
    if yb == zero(T); a = b; end
 
    while b - a > ϵ
        x = (a+b)/2
        y = df(x)
        if y == zero(T)
            a,b = x, x
        elseif sign(y) == sign(ya)
            a = x
        else
            b = x
        end
    end
    return min(a,b)
end

bracket(f, a::Real, b::Real) = bracket(f, promote(a,b)...)
function bracket(f, a::T, b::T) where T <: Real
    fa, fb = f(a), f(b)

    # assuming minimization
    if fa < fb
        b, a = a, b
        fa, fb = fb, fa
    end

    a, fa, b, fb
    c = b + Φ * (b - a)
    fc = f(c)

    while fc < fb

        tmp1 = (b - a) * (fb - fc)
        tmp2 = (b - c) * (fb - fa)

        val = tmp2 - tmp1
        denom = abs(val) < ϵ_min ? 2 * ϵ_min : 2 * val

        w = b - ((b - c) * tmp2 - (b - a) * tmp1) / denom
        wlim = b + growLimit * (c - b)

        local fw
        if (w - c) * (b - w) > zero(T)
            fw = f(w)
            if fw < fc
                a, b = b, w
                fa, fb = fb, fw
                break # IS this needed?
            elseif fw > fb
                c, fc = w, fw
                break # IS this needed?
            end
            w = c + Φ * (c - b)
            fw = f(w)
        elseif (w - wlim) * (wlim - c) >= zero(T)
            w = wlim
            fw = f(w)
        elseif (w - wlim) * (c - w) > zero(T)
            fw = f(w)
            if fw < fc
                b, c = c, w
                w = c + Φ * (c - b)
                fb, fc = fc, fw
                fw = f(w)
            end
        else
            w = c + Φ * (c - b)
            fw = f(w)
        end

        a, fa = b, fb
        b, fb = c, fc
        c, fc = w, fw
    end

    lo, flo   = a, fa
    mid, fmid = b, fb
    hi, fhi   = c, fc

    if lo > hi
        lo, hi = hi, lo
        flo, fhi = fhi, flo
    end

    return (hi => fhi, mid =>fmid, lo => flo)
end