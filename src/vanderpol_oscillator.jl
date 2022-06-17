using DifferentialEquations, CairoMakie, LinearAlgebra, Interpolations, Statistics, DrWatson, NLsolve;

function load_parameters()
    p = Dict()
    p["μ"] = .5
    p["xeq1"] = 0.0
    p["t₁"] = 0.0         # [s] time at which the ramp begins.
    p["F_crit"] = 1                                                 # [N] external force leading to tipping at equilibrium.
    p["F_const"] = p["F_crit"] - 1e-1                               # [N] external force at a point right before tipping for control run.
    return p
end

# 1st-order reformulation derived and double-checked with wikipedia.
function vdP(u, p, t)
    # avoid runaway, as vdP unstable.
    if norm(u, 2) > 10.
        return [0., 0.]
    else
        du1 = u[2]
        du2 = p["μ"] * (1 - u[1]^2) * u[2] - u[1]
        return -[du1, du2]      # we want the formulation with unstable limit cycle, hence the negative sign.
    end
end

function forced_vdP(u, p, t)
    F = get_F(t)
    return vdP(u, p, t) + F .* [0, 1]
end

vdP_stream(u) = Point2f( vdP(u, p, 0.0) )

function solve_forced_vdP(u₀, tspan, p)
    ivp = ODEProblem(forced_vdP, u₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense=true)
    end
end

function solve_vdP(u₀, tspan, p)
    ivp = ODEProblem(vdP, u₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense=true)
    end
end
# fig = Figure(resolution = (800, 800), font = srcdir("cmunrm.ttf"), fontsize = 20)
# ax = Axis(fig[1, 1][1, 1], xlabel = L"$x_{1}$ [m]", ylabel = L"$x_{2}$ [m/s]")
# sp = streamplot!(
#     ax,
#     vdP_stream,
#     -3 .. 3,
#     -3 .. 3,
#     gridsize = (32, 32),
# )

# fig