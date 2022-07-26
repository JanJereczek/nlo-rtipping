include("transfer_func.jl")
include("utils.jl")
mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
end

function get_scatter(node, plot_type, sss, solve_ode)
    Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

    if (plot_type == "Δx") | (plot_type == "single")
        Δx = node
    elseif plot_type == "σ"
        p["σ"] = node
        Δx = 0.6
    elseif plot_type == "D"
        d = D2d(node, p)
        p["d"] = d
        Δx = p["Δx"]
    end

    x₀ = get_x₀(p, Δx)
    println("Getting results for x₀ = $x₀")
    sol_dict = Dict()
    t_buffer = 100

    for a in sss.avec
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
            sol_dict[string(Fmax, "  ", a)] = 0 # sol
            # append!(x_scat, norm(last(sol), 2))
            nend = 100
            end_sol = sol(range(tend - t_buffer / 10, stop = tend, length = nend))
            end_position = reverse([end_sol.u[end-i][1] for i = 0:nend-1])
            append!(x_scat, abs(mean(end_position)))
        end
    end
    return [Fmax_scat, a_scat, x_scat]
end

# function get_scatter(node, plot_type, sss)
#     Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

#     if (plot_type == "Δx") | (plot_type == "single")
#         Δx = node
#     elseif plot_type == "σ"
#         p["σ"] = node
#         Δx = 0.6
#     elseif plot_type == "D"
#         d = D2d(node, p)
#         p["d"] = d
#         Δx = 0.6
#     end

#     x₀ = get_x₀(p, Δx)
#     println("Getting results for x₀ = $x₀")

#     for a in sss.avec
#         p["aF"] = a
#         local tend = p["F_crit"] / a + 10
#         local tspan = (0.0, tend)   # Simulation time span.

#         for Fmax in sss.Fvec

#             p["Fmax"] = Fmax
#             p["t₂"] = p["t₁"] + Fmax / a
#             append!(Fmax_scat, Fmax)
#             append!(a_scat, a)

#             local sol = solve_plo(x₀, tspan, p)
#             append!(x_scat, last(sol)[1])
#         end
#     end
#     return [Fmax_scat, a_scat, x_scat]
# end

function plot_scatter(x, sss, grid_axs, grid_fig, node, prefix_anim, cb)

    if cb == "fixed_cb"
        hm1 = scatter!(
            grid_axs[1],
            x[2],
            x[1],
            color = x[3],
            colormap = sss.cb_maps[1],
            colorrange = sss.cb_limits[1],
        )
        hm2 = scatter!(
            grid_axs[2],
            x[2],
            x[1],
            color = x[4]["1"],              # Only plot γ₁ (or any other γ you wish)
            colormap = sss.cb_maps[2],
            colorrange = sss.cb_limits[2],
        )
    elseif cb == "adaptive_cb"
        hm1 = scatter!(grid_axs[1], x[2], x[1], color = x[4], colormap = sss.cb_maps[1])
        hm2 = scatter!(grid_axs[2], x[2], x[1], color = x[3], colormap = sss.cb_maps[2])
    end

    Colorbar(grid_fig[1, 1][1, 2], hm1)
    Colorbar(grid_fig[1, 2][1, 2], hm2)
    save_fig(prefix_anim, "=$node", "both", grid_fig)
end