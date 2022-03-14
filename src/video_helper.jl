include("transfer_func.jl")
include("utils.jl")
mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
end

function get_scatter(node, plot_type, sss)
    Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

    if (plot_type == "Δx") | (plot_type == "single")
        Δx = node
    elseif plot_type == "σ"
        p["σ"] = node
        Δx = 0.6
    elseif plot_type == "D"
        p["d"] = node
        Δx = 0.6
    end

    x₀ = get_x₀(p, Δx)
    println("Getting results for x₀ = $x₀")

    for a in sss.avec
        p["aF"] = a
        local tend = p["F_crit"] / a + 10
        local tspan = (0.0, tend)   # Simulation time span.

        for Fmax in sss.Fvec

            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + Fmax / a
            append!(Fmax_scat, Fmax)
            append!(a_scat, a)

            local sol = solve_nlo(x₀, tspan, p)
            append!(x_scat, last(sol)[1])
        end
    end
    return [Fmax_scat, a_scat, x_scat]
end

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