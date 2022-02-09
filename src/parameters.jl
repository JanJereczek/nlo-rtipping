include(srcdir("mechanics.jl"))

function load_parameters()

    global p = Dict()
    p["m"] = 10.0         # [kg] mass of oscillator.
    p["c₁"] = 50.0        # [kg/s²] stiffness of linear spring.
    p["d"] = 20.0         # [kg/s] dampening factor. Criteria for non-oscillating behavior: d > 2*sqrt(m*c)
    p["k₁"] = 50.0        # [kg/s²] max stiffness of non-linear spring.
    p["k₂"] = 1e3         # [1] steepness param of sigmoid.
    p["xₜ"] = 1.5         # [m] tipping displacement.
    p["t₁"] = 0.0         # [s] time at which the ramp begins.
    p["g"] = 10.0         # [m/s²] earth acceleration constant. (simplified number to have round thresholds).

    p["F_crit"] = (p["c₁"] + p["k₁"]) * p["xₜ"] - p["m"] * p["g"]   # [N] external force leading to tipping at equilibrium.
    p["F_const"] = p["F_crit"] - 1e-1                               # [N] external force at a point right before tipping for control run.
    p["F1"] = (p["c₁"] + p["k₁"]) * p["xₜ"]                         # [N] total force leading to tipping at equilibrium.

    p["D"] = get_D(p)
    p["xeq1"] = p["m"] * p["g"] / (p["c₁"] + p["k₁"])
    p["xeq2"] = p["m"] * p["g"] / p["c₁"]
    return p
end
