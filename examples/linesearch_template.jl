"""
 Algorithm explanation
"""


 @with_kw struct SampleLineSearch{T}
    f_tol::T = 1e-4 # c_1 Wolfe sufficient decrease condition
    gtol::T = 0.9   # c_2 Wolfe curvature condition (Recommend 0.1 for GradientDescent)
    x_tol::T = 1e-8
    alphamin::T = 1e-16
    alphamax::T = 65536.0
    maxfev::Int = 100
end

 
 function (ls::SampleLineSearch)(df::AbstractObjective, x::AbstractArray{T},
                           s::AbstractArray{T}, α::Real,
                           x_new::AbstractArray{T}, phi_0::Real, dphi_0::Real) where T
     ϕ, ϕdϕ = make_ϕ_ϕdϕ(df, x_new, x, s)
     # or
     # ϕdϕ = make_ϕdϕ(df, x_new, x, s)

     ls(ϕdϕ, alpha, ϕ_0, dϕ_0)
    end
 
 (ls::SampleLineSearch)(ϕ, dϕ, ϕdϕ, phi_0, dphi_0) = ls(ϕ, ϕdϕ, phi_0, dphi_0)
 
 # TODO: Should we deprecate the interface that only uses the ϕ and ϕd\phi arguments?
 function (ls::SampleLineSearch)(ϕ, ϕdϕ,
                           phi_0::Real,
                           dphi_0::Real) where T # Should c and phi_0 be same type?
            @unpack f_tol, gtol, x_tol, alphamin, alphamax, maxfev = ls
 
    # Do stuff ... 

    # Include optimizer for argmin_α f(x + α * sn)

    return  α, sn # It must return the step size α and the search direction at step_n
 end