include("transfer_func.jl")

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

    x₁_t₀ = p["xeq1"] - Δx
    x₀ = [x₁_t₀, 0.0]
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

            Fstep = -get_c(x₁_t₀) * Δx
            ft_F = ft_stepramp(p["t₁"], p["t₂"], Fmax, Fstep, a, sss.ω_vec)

            get_ft_F = LinearInterpolation(sss.ω_vec, abs.(ft_F))
            # append!(resfreq_scat, get_ft_F(ω_res[2]) )
            # ωres_sample = range(sss.ω_res[1], stop = sss.ω_res[3], length = sss.n_int)
            # ωreq_product = get_ft_F.(ωres_sample) .* get_ω_ampresp.(ωres_sample)
            # A_ω_res_integrated = sum(ωreq_product) / sss.n_int

            # ωreq_product = get_ft_F.(sss.ω_vec) .* get_ω_ampresp.(sss.ω_vec)
            # A_ω_res_integrated = sum(ωreq_product) / length(sss.ω_vec)


            #tf = nlo_transfer1(p, x₀, sss.ω_vec)
            #append!(resfreq_scat, maximum(abs.(tf)))

            # fig = Figure(resolution = (1000, 500))
            # ax = Axis(
            #     fig[1, 1],
            #     xlabel = "ω [rad/s]",
            #     ylabel = "|Y(jw)|",
            #     xscale = log10,
            #     yscale = log10,
            #     xminorticks = IntervalsBetween(9),
            #     xminorticksvisible = true,
            #     xminorgridvisible = true,
            # )
            # lines!(ax, sss.ω_vec, abs.(tf))
            # xlims!(ax, (1e-3, 1e3))
            # ylims!(ax, (1e-3, 1e3))
            # save_fig("plots/", "hash=$(randn(1)[1])", "png", fig)

            # tf_res = nlo_transfer1(p, x₀, sss.ω_res[2] )
            # append!(resfreq_scat, abs.(tf_res))

            tf_res = nlo_transfer1(p, x₀, range(sss.ω_res[1], stop = sss.ω_res[3], length=1000) )
            append!(resfreq_scat, sum(abs.(tf_res)))
            # append!(resfreq_scat, sum(abs.(log10.(tf_res))))
        end
    end
    return [Fmax_scat, a_scat, resfreq_scat, x_scat]
end

function plot_scatter_fixedcb(x, sss, grid_axs, grid_fig, node, prefix_anim)
    # println(sum(x[3]), sum(x[4]))
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
    hm1 = scatter!(grid_axs[1], x[2], x[1], color = x[4], colormap = sss.cb_maps[1])

    hm2 = scatter!(grid_axs[2], x[2], x[1], color = x[3], colormap = sss.cb_maps[2])
    Colorbar(grid_fig[1, 1][1, 2], hm1)
    Colorbar(grid_fig[1, 2][1, 2], hm2)
    save_fig(prefix_anim, "=$node", "both", grid_fig)
end
