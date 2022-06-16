include("mechanics.jl")
include("forcing.jl")

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

function nl_osc_free(x)
    c = get_c(x[1])
    return [x[2], -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + p["g"]]
end

nl_osc_free_stream(x) = Point2f( nl_osc_free(x) )
nl_osc_free_quiver(x) = 5e-2 .* nl_osc_free(x)

function nl_oscillator!(dx, x, p, t)
    F = get_F(t)
    dx = nl_osc_free(x) + F .* [0, 1/( p["m"] )]
    # dx[1] = x[2]
    # dx[2] = -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + F / p["m"] + p["g"]
end

function nl_oscillator_F!(dx, x, p, t)
    F = get_F(t)
    dx = nl_osc_free(x) .+ F .* [0, 1/( p["m"] )] .- [0, p["g"]]
    # dx[1] = x[2]
    # dx[2] = -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + F / p["m"]
end

function solve_nlo(x₀, tspan, p)
    ivp = ODEProblem(nl_oscillator!, x₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense=true)
    end
end

function solve_nlo_F(x₀, tspan, p)
    ivp = ODEProblem(nl_oscillator_F!, x₀, tspan, p)
    if p["fixed_dt"]
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp, dense=true)
    end
end