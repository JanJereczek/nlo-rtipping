include("utils.jl")

# Transfer function of the system in the left half of the state space
function nlo_transfer(p::Dict, x₀::Vector, ω::Any, region::Int)
    s = im .* ω
    if region == 1
        c = p["c₁"] + p["k₁"]
    elseif region == 2
        c = p["c₁"]
    end

    K = 1 / p["m"]
    denum = s .^ 2 .+ 2*p["D"]*p["ω₀1"] .* s .+ p["ω₀1"]^2
    G₁ = K ./ denum

    num_ic = x₀[1] .* s .+ x₀[2] .+ 2*p["D"]*p["ω₀1"] * x₀[1]
    G₂ = num_ic ./ denum
    return G₁, G₂
end


function get_Y(G₁, G₂, U)
    return G₁ .* U .+ G₂
end