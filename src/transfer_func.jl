include("utils.jl")

function nlo_transfer1(p, x₀, ω)
    s = im .* ω
    c = p["c₁"] + p["k₁"]
    denum = p["m"] .* s .^ 2 .- p["d"] .* s .- c
    num_ic = (p["m"] .* s .+ p["d"]) .* x₀[1] .+ p["d"] * x₀[2]
    u = ft_stepramp(p["t₁"], p["t₂"], p["Fmax"], p["m"]*p["g"], p["aF"], ω)
    tf = num_ic ./ denum .+ 1 ./ denum .* u
    return tf
end