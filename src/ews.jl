
function get_bassin_boundary(p, Fconst)
    # Load the important constants
    w0 = p["ω₀1"]
    wd = p["ωd"]
    D = p["D"]
    xₜ = p["xₜ"]

    # Define the constant forcing
    p["F_type"] = "const"
    p["F_const"] = Fconst
    p["fixed_dt"] = true
    p["dt"] = 0.01
    xmg = (p["m"] * p["g"] + p["F_const"]) / (p["c₁"] + p["k₁"])

    # Analytical solutions for position and speed.
    x₁(A, B, t) = exp(-D*w0*t) * ( A*cos(wd*t) + B*sin(wd*t) ) + xmg
    x₂(A, B, t) = -D*w0*exp(-D*w0*t) * ( A*cos(wd*t) + B*sin(wd*t) ) + exp(-D*w0*t) * ( -A*wd*sin(wd*t) + B*wd*cos(wd*t) )
    
    # Analytical solutions for the constants A and B. The latter depends on an unknown time t_tip.
    # Equations derived from x₁(t = 0) = x₁(t = t_tip) = xₜ
    A = xₜ - xmg
    fB(t) = (exp(D*w0*t) * (xₜ - xmg) - A*cos(wd*t)) / (sin(wd*t))

    # Compute t_tip by setting x₂(t = t_tip) = 0. Notice: Initial guess quite important... but not dependent on the input Fconst!
    function f!(F, x)
        F[1] = x₂(A, fB(x[1]), x[1])
    end
    lb, ub, initial_x = [1.0], [5.0], [2.0]
    @time sol_nl = mcpsolve(f!, lb, ub, initial_x, reformulation = :smooth)
    B = fB(t_tip)

    # Get the trajectory on bassin boundary.
    x₀ = [x₁(A, B, 0), x₂(A, B, 0)]
    return solve_nlo(x₀, (0, t_tip), p)
end