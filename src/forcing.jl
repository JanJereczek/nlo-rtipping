function get_F(t)

    if p["F_type"] == "const"
        F = p["F_const"]
    elseif t < p["t₁"]
        F = 0
    elseif t < p["t₂"]
        F = p["aF"] * (t - p["t₁"])
    else
        F = p["aF"] * (p["t₂"] - p["t₁"])
    end

    if p["F_type"] == "noisy"
        F += p["σ"] * randn(1)[1]
    end
    return F
end

function get_Fvec(t)
    F = similar(t)
    F[t.<p["t₁"]] .= 0
    indices = (t .< p["t₂"]) .& (t .≥ p["t₁"])
    F[indices] .= p["aF"] * (t[indices] .- p["t₁"])
    F[t.≥p["t₂"]] .= p["aF"] * (p["t₂"] - p["t₁"])
    return F
end
