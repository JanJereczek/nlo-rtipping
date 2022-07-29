include("transfer_func.jl")
include("utils.jl")
mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
end

function get_scatter!(node, plot_type, sss, solve_ode, grid_dict, sol_dict)
    println( string("Getting results for ", plot_type, " = $node") )
    Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

    if (plot_type == "Δx")
        Δx = node
    elseif plot_type == "D"
        d = D2d(node, p)
        p["d"] = d
        Δx = p["Δx"]
    end

    x₀ = get_x₀(p, Δx)
    t_buffer = 100

    @showprogress for a in sss.avec
        p["aF"] = a
        if p["fixed_tend"] == true
            local tend = p["F_crit"] / minimum(sss.avec) + t_buffer
        else
            local tend = p["F_crit"] / a + t_buffer
        end
        local tspan = (0.0, tend)   # Simulation time span.

        for Fmax in sss.Fvec

            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + Fmax / a
            append!(Fmax_scat, Fmax)
            append!(a_scat, a)

            local sol = solve_ode(x₀, tspan, p)
            # append!(x_scat, norm(last(sol), 2))
            # sol_dict[string(Fmax, "  ", a)] = sol
            nend = 100
            end_sol = sol(range(tend - t_buffer / 10, stop = tend, length = nend))
            end_position = reverse([end_sol.u[end-i][1] for i = 0:nend-1])
            append!(x_scat, abs(mean(end_position)))
        end
    end
    grid_dict[string(node)] = [Fmax_scat, a_scat, x_scat]
end
