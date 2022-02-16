include("transfer_func.jl")
include("utils.jl")
mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
    ω_vec::Vector{Float64}
    ω_res::Vector{Float64}
    n_int::Int64
    cb_limits::Vector{Tuple{Float64,Real}}
    cb_maps::Vector{Symbol}
end

function get_scatter(node, anim_type, sss)
    Fmax_scat, a_scat, resfreq_scat, x_scat = zeros(0), zeros(0), zeros(0), zeros(0)

    if (anim_type == "Δx") | (anim_type == "none")
        Δx = node
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

            Fstep = get_c(x₀[1]) * Δx
            ft_F = ft_stepramp(p["t₁"], p["t₂"], Fstep, a, sss.ω_vec)

            get_ft_F = LinearInterpolation(sss.ω_vec, abs.(ft_F))
            ω_res = range(sss.ω_res[1], stop = sss.ω_res[3], length=1000)
            G₁_res, G₂_res = nlo_transfer(p, x₀, ω_res, 1)
            U_res = get_ft_F(ω_res)
            Y_res = get_Y(G₁_res, G₂_res, U_res)
            append!(resfreq_scat, sum(abs.(Y_res)))
            # append!(resfreq_scat, sum(abs.(log10.(tf_res))))
        end
    end
    return [Fmax_scat, a_scat, resfreq_scat, x_scat]
end

function plot_scatter_fixedcb(x, sss, grid_axs, grid_fig, node, prefix_anim)
    # x[3] = normalise(x[3])
    # println( extrema(x[3]) )

    scatter!(
        grid_axs[1],
        x[2],
        x[1],
        color = x[4],
        colormap = sss.cb_maps[1],
        colorrange = sss.cb_limits[1],
    )
    scatter!(
        grid_axs[2],
        x[2],
        x[1],
        color = x[3],
        colormap = sss.cb_maps[2],
        colorrange = sss.cb_limits[2],
    )
    save_fig(prefix_anim, "=$node", "both", grid_fig)
end

function plot_scatter(x, sss, grid_axs, grid_fig, node, prefix_anim)
    # println(sum(x[3]), sum(x[4]))
    dist2mean = normalise(x[3])

    hm1 = scatter!(grid_axs[1], x[2], x[1], color = x[4], colormap = sss.cb_maps[1])

    hm2 = scatter!(grid_axs[2], x[2], x[1], color = dist2mean, colormap = sss.cb_maps[2])
    Colorbar(grid_fig[1, 1][1, 2], hm1)
    Colorbar(grid_fig[1, 2][1, 2], hm2)
    save_fig(prefix_anim, "=$node", "both", grid_fig)
end
