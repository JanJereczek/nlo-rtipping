# Define parameters dictionary with standard values.
function standard_parameters()
    global p = Dict()
    p["m"] = 10.0         # [kg] mass of oscillator.
    p["c₁"] = 50.0        # [kg/s²] stiffness of linear spring.
    p["d"] = 20.0         # [kg/s] dampening factor. Criteria for non-oscillating behavior: d > 2*sqrt(m*c)
    p["xₜ"] = 1.5         # [m] tipping displacement.
    return p
end

# Compute resonance frequency (in rad/s)
get_ω₀(c) = sqrt(c / p["m"])

# Compute damping degree D.
get_D(p, x::Float64) = p["d"] / (2 * sqrt(p["m"] * get_total_plastic_stiffness(x)))

# Compute sigmoid with offset xi (inflection point) and steepness scaling α
parametric_sigmoid(α::Float64, xi::Float64, x::Float64) = 1. - 1. / (1. + exp(- α * (x - xi)) )

# Compute parametric sigmoid for given values of xi and α
sigmoid(x::Float64) = parametric_sigmoid(α2, x_inflection_spring2, x)

# Compute stiffness of plastic spring 2.
get_plastic_stiffness(x::Float64) = c2_min + Δc2 * sigmoid(x)

# Compute total stiffness of elastic spring 1 and plastic spring 2.
get_total_plastic_stiffness(x::Float64) = p["c₁"] .+ get_plastic_stiffness(x)

# Compute force needed to maintain system at equilibrium x.
po_equilibrium_F(x::Float64) = get_total_plastic_stiffness(x) * x - p["m"] * g

function plastic_oscillator(x::Vector{Float64})
    c = get_total_plastic_stiffness(x[1])
    return [x[2], -c / p["m"] * x[1] - p["d"] / p["m"] * x[2] + g]
end

function forced_plastic_oscillator(x, p, t)
    F = p["forcing"](t)
    return plastic_oscillator(x) + F .* [0, 1 / (p["m"])]
end

plastic_oscillator_stream(x) = Point2f(plastic_oscillator(x))
plastic_oscillator_quiver(x) = 5e-2 .* plastic_oscillator(x)
solve_po(x₀, tspan, p) = solve(ODEProblem(forced_plastic_oscillator, x₀, tspan, p), dense = true)

# Compute Bode diagram of po.
function plastic_amp_response(ω::Float64, x::Float64)
    c = get_total_plastic_stiffness(x)
    ω₀ = get_ω₀( c )
    η = ω ./ ω₀
    D = get_D(p, x)
    return 1 / (c * (sqrt((1 - η^2)^2 + 4 * D^2 * η^2)))
end

# Test
# x1 = 0:.01:2
# lines( x1, get_plastic_stiffness.(x1) )
# lines( x1, cumsum(get_plastic_stiffness.(x1)) )