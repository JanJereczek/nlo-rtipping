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
    G = K ./ denum

    num_ic = x₀[1] .* s .+ x₀[2] .+ 2*p["D"]*p["ω₀1"] * x₀[1]
    Y₀ = num_ic ./ denum
    return G, Y₀
end


function get_Y(G, Y₀, U)
    return G .* U .+ Y₀
end

function get_Y_stat(x₀, ω)
    return x₀ ./ (im * ω)
end