include(srcdir("parameters.jl"))

function get_analytical_sol_IC(t)
    k1 = exp( -p["D"] * t * p["ω₀1"] )
    k2 = p["D"] * p["ω₀1"] * p["x0"] + p["dx0"]
    k3 = sqrt( p["D"]^2 -1 ) * p["ω₀1"]
    return k1 * ( k2 * sinh( k3*t ) / k3 + p["x0"] * cosh( k3*t )  )
end

p = load_parameters()
p["x0"] = 0.4
p["dx0"] = 0