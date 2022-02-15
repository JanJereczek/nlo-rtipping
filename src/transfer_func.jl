include("utils.jl")

# Transfer function of the system in the left half of the state space
function nlo_transfer1(p, x₀, ω)
    s = im .* ω
    c = p["c₁"] + p["k₁"]

    K = 1 / p["m"]
    denum = s .^ 2 .+ 2*p["D"]*p["ω₀1"] .* s .+ p["ω₀1"]^2
    G₁ = K ./ denum

    num_ic = x₀[1] .* s .+ x₀[2] .+ 2*p["D"]*p["ω₀1"] * x₀[1]
    G₂ = num_ic ./ denum

    U = ft_stepramp(p["t₁"], p["t₂"], p["Fmax"], p["m"]*p["g"], p["aF"], ω)
    # U = ft_stepramp(p["t₁"], p["t₂"], p["Fmax"], 0.0, p["aF"], ω)
    Y = G₁ .* U .+ G₂

    return G₁, G₂, U, Y
end