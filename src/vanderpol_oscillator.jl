using DifferentialEquations,
    CairoMakie, LinearAlgebra, Interpolations, Statistics, DrWatson, NLsolve;

function load_parameters()
    p = Dict()
    p["μ"] = 0.5
    p["xeq1"] = 0.0
    p["t₁"] = 0.0           # [s] time at which the ramp begins.
    p["F_crit"] = 1         # [N] results from the bifurcation analysis.
    p["F_const"] = p["F_crit"] - 1e-1                               # [N] external force at a point right before tipping for control run.
    return p
end

# 1st-order reformulation derived and double-checked with wikipedia.
function vdP(u, p, t)
    # unstable limit cycle --> negative sign.
    du1 = -( u[2] )
    du2 = -( p["μ"] * (1 - u[1]^2) * u[2] - u[1] )
    return [du1, du2]
end

function forced_vdP(u, p, t)
    # avoid runaway, as vdP unstable.
    if norm(u, 2) > 10.0
        return [0.0, 0.0]
    else
        F = get_F(t)
        return vdP(u, p, t) - F .* [0, 1]
    end
end

constant_forced_vdP(u, p, t) = vdP(u, p, t) - p["F₀"] .* [0, 1]
vdP_stream(u) = Point2f(vdP(u, p, 0.0))
forced_vdP_stream(u, F) = Point2f(vdP(u, p, 0.0) - F .* [0, 1])

function solve_vdP(u₀, tspan, p)
    ivp = ODEProblem(vdP, u₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense = true)
    end
end

@inline function solve_forced_vdP(u₀, tspan, p)
    ivp = ODEProblem(forced_vdP, u₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense = true)
    end
end

function solve_constantly_forced_vdP(u₀, tspan, p)
    ivp = ODEProblem(constant_forced_vdP, u₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense = true)
    end
end
