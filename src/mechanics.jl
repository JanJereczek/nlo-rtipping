function get_nl_stiffness(x₁)
    return p["k₁"] / (1 + exp(p["k₂"] * (x₁ - p["xₜ"])))
end

function get_D(p)
    return p["d"] / (2 * sqrt(p["m"] * (p["c₁"] + p["k₁"])))
end

function get_bode_amp(η)
    D = get_D(p)
    return (sqrt((1 - η^2)^2 + 4 * D^2 * η^2))^(-1)
end

function get_resonance_freq(p)
    ω_res = sqrt(1 - 2 * get_D(p) .^ 2)
    f_res = ω_res / (2 * π)
    return f_res, ω_res
end
