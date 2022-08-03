include("forcing.jl")
include("mechanics_plo.jl")

function standard_parameters()
    global p = Dict()
    p["m"] = 10.0         # [kg] mass of oscillator.
    p["c₁"] = 50.0        # [kg/s²] stiffness of linear spring.
    p["d"] = 20.0         # [kg/s] dampening factor. Criteria for non-oscillating behavior: d > 2*sqrt(m*c)
    p["k₁"] = 50.0        # [kg/s²] max stiffness of non-linear spring.
    p["k₂"] = 1e6         # [1] steepness param of sigmoid.
    p["xₜ"] = 1.5         # [m] tipping displacement.
    p["t₁"] = 0.0         # [s] time at which the ramp begins.
    p["g"] = 10.0         # [m/s²] earth acceleration constant. (simplified number to have round thresholds).
    return p
end

function compute_parameters(p)
    p["ω₀1"] = get_ω₀(p["c₁"] + p["k₁"])
    p["ω₀2"] = get_ω₀(p["c₁"])
    p["D"] = get_D(p)
    p["ωd"] = p["ω₀1"] * sqrt(1 - p["D"]^2)
    p["xeq1"] = p["m"] * p["g"] / (p["c₁"] + p["k₁"])
    p["xeq2"] = p["m"] * p["g"] / p["c₁"]

    p["F_crit"] = (p["c₁"] + p["k₁"]) * p["xₜ"] - p["m"] * p["g"]   # [N] external force leading to tipping at equilibrium.
    p["F_const"] = p["F_crit"] - 1e-1                               # [N] external force at a point right before tipping for control run.
    p["F1"] = (p["c₁"] + p["k₁"]) * p["xₜ"]                         # [N] total force leading to tipping at equilibrium.
    return p
end

function load_parameters()
    p = standard_parameters()
    p = compute_parameters(p)
    return p
end

function load_highfreq_parameters()

    p = load_parameters()
    p["m"] = 1.0            # [kg] mass of oscillator.
    # p["xₜ"] = 0.15           # [m] tipping displacement.
    p["d"] = 1.0            # [kg/s] dampening factor. Criteria for non-oscillating behavior: d > 2*sqrt(m*c)
    p = compute_parameters(p)
    return p
end

function branch_xeq1(F)
    if F < p["F1"]
        return F / (p["c₁"] + p["k₁"])
    else
        return F / (p["c₁"])
    end
end

function branch_xeq2(F)
    return F / (p["c₁"])
end

function get_c(x1)
    return p["c₁"] .+ get_nl_stiffness.(x1)
end

function pwlin_oscillator(x)
    c = get_c(x[1])
    return [x[2], -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + p["g"]]
end

pwlin_oscillator_stream(x) = Point2f(pwlin_oscillator(x))
pwlin_oscillator_quiver(x) = 5e-2 .* pwlin_oscillator(x)

function forced_pwlin_oscillator(x, p, t)
    F = get_F(t)
    return pwlin_oscillator(x) + F .* [0, 1 / (p["m"])]
end

function nogravity_forced_pwlin_oscillator(x, p, t)
    F = get_F(t)
    return pwlin_oscillator(x) .+ F .* [0, 1 / (p["m"])] .- [0, p["g"]]
end

function solve_plo(x₀, tspan, p)
    ivp = ODEProblem(forced_pwlin_oscillator, x₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, BS3(), dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense = true)
    end
end

function solve_plo_F(x₀, tspan, p)
    ivp = ODEProblem(nogravity_forced_pwlin_oscillator, x₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense = true)
    end
end
