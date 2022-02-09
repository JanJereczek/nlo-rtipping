using Interpolations

# Nonlinear stiffness c₂(x).
function get_nl_stiffness(x₁)
    return p["k₁"] / (1 + exp(p["k₂"] * (x₁ - p["xₜ"])))
end

# Compute damping degree D.
function get_D(p)
    return p["d"] / (2 * sqrt(p["m"] * (p["c₁"] + p["k₁"])))
end

# Compute damping factor of damping element based on damping degree and stiffnesses.
function D2d(D, p)
    return D .* (2 * sqrt(p["m"] * (p["c₁"] + p["k₁"])))
end

# Compute Bode diagram of nlo.
function get_bode_amp(η)
    D = get_D(p)
    return (sqrt((1 - η^2)^2 + 4 * D^2 * η^2))^(-1)
end

# Compute resonance frequency for given parameter set.
function get_resonance_freq(p)
    ω_res = sqrt(1 - 2 * get_D(p) .^ 2)
    return ω_res
end

# Compute characteristics of resonance
function get_resonance_characteristics(p::Dict, ω_vec::Vector, res_threshold::Float64)
    amp_resp = get_bode_amp.(ω_vec)
    get_ω_ampresp = LinearInterpolation(ω_vec, amp_resp)        # Create linear interpolations
    get_f_ampresp = LinearInterpolation(ω2f(ω_vec), amp_resp)

    ω_res = round.(get_resonance_freq(p); digits = 3)           # Obtained by deriving bode
    ix_res = argmax(amp_resp)                                   # Index corresponding to resonance
    ω_res_num = ω_vec[ix_res]                                   # Corresponding numerical estimation
    dω = ω_vec[ix_res+1] - ω_vec[ix_res]                        # Frequency stepping
    test_msg = "Analytical and numerical computation of resonance frequency coincide? "
    println(test_msg, isapprox(ω_res, ω_res_num; atol = 1 * dω))

    # Get lower and upper limits of resonance region
    ωres_vec = ω_vec[amp_resp.>res_threshold]
    ω1_res, ω2_res = round(ωres_vec[1]; digits = 3), round(ωres_vec[end]; digits = 3)

    # Zip characteristics
    ω_res_char = [ω1_res, ω_res, ω2_res]
    f_res_char = ω2f(ω_res_char)
    bound_msg = " - lower bound, peak and upper bound: "
    println(string("Resonance freq (rad/s)", bound_msg), ω1_res, "  ", ω_res, "  ", ω2_res)
    println(
        string("Resonance freq (Hz)", bound_msg),
        f_res_char[1],
        "  ",
        f_res_char[2],
        "  ",
        f_res_char[3],
    )

    return get_ω_ampresp, get_f_ampresp, amp_resp, ω_res_char, f_res_char
end
