include("utils.jl")

# Transfer function of the system in the left half of the state space
function nlo_transfer(p::Dict, x₀::Vector, ω::Any, region::Int)
    s = im .* ω
    if region == 1
        c = p["c₁"] + p["k₁"]
        xeq = p["xeq1"]
    elseif region == 2
        c = p["c₁"]
        xeq = p["xeq2"]
    end

    K = 1 / p["m"]
    denum = s .^ 2 .+ 2 * p["D"] * p["ω₀1"] .* s .+ p["ω₀1"]^2
    G = K ./ denum

    # num_ic = x₀[1] .* s .+ x₀[2] .+ 2*p["D"]*p["ω₀1"] * x₀[1]
    num_ic = (x₀[1] - xeq) .* s .+ x₀[2] .+ 2 * p["D"] * p["ω₀1"] * (x₀[1] - xeq)
    Y₀ = num_ic ./ denum
    return G, Y₀
end

# Transfer function of the system in the left half of the state space
function fsplit_transfer(p::Dict, x₀::Vector, ω::Any, U::Vector, region::Int)
    s = im .* ω
    if region == 1
        c = p["c₁"] + p["k₁"]
        xeq = p["xeq1"]
    elseif region == 2
        c = p["c₁"]
        xeq = p["xeq2"]
    end

    K = 1 / p["m"]
    denum = s .^ 2 .+ 2 * p["D"] * p["ω₀1"] .* s .+ p["ω₀1"]^2
    G = K ./ denum
    Y₁ = G .* U

    num_ic = (x₀[1] - xeq) .* s .+ x₀[2] .+ 2 * p["D"] * p["ω₀1"] * (x₀[1] - xeq)
    Y₀ = num_ic ./ denum

    Y₂ = p["g"] ./ s ./ denum
    return G, Y₀, Y₁, Y₂
end

function get_Y(G, Y₀, U)
    return G .* U .+ Y₀
end

function get_Y_stat(x₀, ω)
    return x₀ ./ (im * ω)
end