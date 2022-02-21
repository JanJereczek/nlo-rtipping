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
    Fmax_scat, a_scat, resfreq_scat, x_scat = zeros(0), zeros(0), zeros(0)
    γ₁_scat, γ₂_scat, γ₃_scat, γ₄_scat, γ₅_scat = zeros(0), zeros(0), zeros(0), zeros(0), zeros(0)

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

            # Fstep = get_c(x₀[1]) * Δx
            U = ft_stepramp(p["t₁"], p["t₂"], p["m"]*p["g"], a, sss.ω_vec)
            G, Y₀ = nlo_transfer(p, x₀, sss.ω_vec, 1)
            Y = get_Y(G, Y₀, U)
            Yabs = abs.(Y)
            Ystatabs = abs.(p["Ystat"])

            ix_ω_res = get_box_indices(sss.ω_vec, sss.ω_res[1], sss.ω_res[3])
            ix_ω_t₂ = get_box_indices(sss.ω_vec, π/p["t₂"], 4*π/p["t₂"])
            Y_res = Y[ix_ω_res]
            Y_t₂ = Y[ix_ω_t₂]

            # get_U = LinearInterpolation(sss.ω_vec, abs.(U))
            # ω_res = range(sss.ω_res[1], stop = sss.ω_res[3], length=1000)
            # G_res, Y₀_res = nlo_transfer(p, x₀, ω_res, 1)
            # U_res = get_U(ω_res)
            # Y_res = get_Y(G_res, Y₀_res, U_res)

            # ω_t₂ = range(π/p["t₂"], 4*π/p["t₂"], length=1000)
            # G_t₂, Y₀_t₂ = nlo_transfer(p, x₀, ω_t₂, 1)
            # U_t₂ = get_U(ω_t₂)
            # Y_t₂ = get_Y(G_t₂, Y₀_t₂, U_t₂)
            # γ₁ = sum( U_res )

            ΔY = ( Yabs - Ystatabs ) ./ Ystatabs
            γ₁ = sum( abs.(Y_res) )
            γ₂ = sum( ΔY[ΔY .> 0] )
            γ₃ = maximum(ΔY)
            γ₄ = 1/(2*π) * sum( Yabs.^2 - Ystatabs.^2 )
            γ₅ = 1/(2*π) * sum( abs.(Y_t₂).^2 - Ystatabs.^2 )

            append!(γ₁_scat, γ₁)
            append!(γ₂_scat, γ₂)
            append!(γ₃_scat, γ₃)
            append!(γ₄_scat, γ₄)
            append!(γ₅_scat, γ₅)
        end
    end
    return [Fmax_scat, a_scat, resfreq_scat, x_scat]
end

function plot_scatter(x, sss, grid_axs, grid_fig, node, prefix_anim, cb)
    # println(sum(x[3]), sum(x[4]))
    # dist2mean = normalise(x[3])
    # x[3]

    if cb == "fixed_cb"
        hm1 = scatter!(
            grid_axs[1],
            x[2],
            x[1],
            color = x[4],
            colormap = sss.cb_maps[1],
            colorrange = sss.cb_limits[1],
        )
        hm2 = scatter!(
            grid_axs[2],
            x[2],
            x[1],
            color = x[3],
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
