include("mechanics.jl")
include("forcing.jl")

function xeq1(F)
    if F < p["F2"]
        return F / (p["c₁"] + p["k₁"])
    else
        return F / (p["c₁"])
    end
end

function xeq2(F)
    if F < p["F1"]
        return F / (p["c₁"] + p["k₁"])
    else
        return F / (p["c₁"])
    end
end

function get_c(x1)
    return p["c₁"] .+ get_nl_stiffness.(x1)
end

function nl_oscillator!(dx, x, p, t)
    c = get_c(x[1])
    F = get_F(t)
    dx[1] = x[2]
    dx[2] = -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + F / p["m"] + p["g"]
end

function nl_osc_free_quiver(x1, x2)
    scale = 5e-2
    c = get_c(x[1])
    return [scale .* x2, scale .* (-c ./ p["m"] .* x1 - p["d"] / p["m"] .* x2 .+ p["g"])]
end

function nl_osc_free(x)
    c = get_c(x[1])
    return [x[2], -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + p["g"]]
end

function nl_osc_free_stream(x1, x2)
    c = get_c(x1)
    return Point(x2, -c / p["m"] * x1 - p["d"] / p["m"] * x2 + p["g"])
end

# function get_sol(v₀, tspan, p)
#   local x₀ = [0.0, v₀]      # Set an initial speed for the system to move towards tipping.
#   local ivp = ODEProblem(nl_oscillator!, x₀, tspan, p)
#   return solve(ivp)
# end

function solve_nlo(x₀, tspan, p)
    # println(x₀)
    ivp = ODEProblem(nl_oscillator!, x₀, tspan, p)
    if p["F_type"] == "noisy"
        return solve(ivp, dt = p["dt"], adaptive = false)
    else
        return solve(ivp)
    end
end
