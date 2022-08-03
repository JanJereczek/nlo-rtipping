
function get_bode_a(a, p, Δx, Fvec, ω_vec)
    tf = Dict()
    tf["x₀"] = get_x₀(p, Δx)
    p_temp = p
    j = 28
    j = 32
    p_temp["Fmax"], p_temp["aF"], p_temp["t₂"] = Fvec[j], a, Fvec[j] / a

    tf["U"] = ft_stepramp(
        p_temp["t₁"],
        p_temp["t₂"],
        p_temp["m"] * p_temp["g"],
        p_temp["aF"],
        ω_vec,
    )
    tf["G"], tf["Y₀"] = nlo_transfer(p_temp, tf["x₀"], ω_vec, region)
    tf["Y"] = get_Y(tf["G"], tf["Y₀"], tf["U"])
    tf["Ystat"] = p["Ystat"]
    println("a=", round(a; digits = 2), ":  ", sum(abs.(tf["Y"]) .> abs.(tf["Ystat"])))
    plot_bode(p, ω_vec, string(prefix, "a=", a), tf)
end
