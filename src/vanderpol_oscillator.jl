using DifferentialEquations, CairoMakie, Colors, Interpolations, Statistics, DrWatson, NLsolve;

function load_parameters()
    p = Dict()
    p["μ"] = 0.5
    p["F_crit"] = (p["c₁"] + p["k₁"]) * p["xₜ"] - p["m"] * p["g"]   # [N] external force leading to tipping at equilibrium.
    p["F_const"] = p["F_crit"] - 1e-1                               # [N] external force at a point right before tipping for control run.
    p["F1"] = (p["c₁"] + p["k₁"]) * p["xₜ"]                         # [N] total force leading to tipping at equilibrium.
    return p
end

# 1st-order reformulation derived and double-checked with wikipedia.
function vdP(x, p, t)
    dx1 = x[2]
    dx2 = p["μ"] * (1 - x[1]^2) * x[2] - x[1]
    return -[dx1, dx2]      # we want the formulation with unstable limit cycle, hence the negative sign.
end

function vdP!(dx, x, p, t)
    F = get_F(t)
    dx = vdp(x, p, t) + F .* [0, 1]
end

vdP_stream(x) = Point2f( vdP(x, p, 0.0) )

function solve_vdP(x₀, tspan, p)
    ivp = ODEProblem(vdP!, x₀, tspan, p)
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