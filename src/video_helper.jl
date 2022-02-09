mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
    ω_vec::Vector{Float64}
    f_vec::Vector{Float64}
    f1_res::Float64
    f2_res::Float64
    n_int::Int64
    cb_limits::Vector{Tuple{Float64, Real}}
    cb_maps::Vector{Symbol}
end

function get_scatter(Δx, sss)
    Fmax_scat, a_scat, resfreq_scat, x_scat = zeros(0), zeros(0), zeros(0), zeros(0)
    x0 = p["xeq1"] - Δx
    x₀ = [x0, 0.0]
    println("Getting results for x₀ = $x₀")

    for a in sss.avec
        p["aF"] = a
        local tend = p["F_crit"] / a + 10
        local tspan = (0.0, tend)   # Simulation time span.

        for Fmax in sss.Fvec

            p["F_max"] = Fmax
            p["t₂"] = p["t₁"] + Fmax / a
            append!(Fmax_scat, Fmax)
            append!(a_scat, a)

            local sol = solve_nlo(x₀, tspan, p)
            append!(x_scat, last(sol)[1])
            # println(x_scat)

            # ft_F = ft_ramp(p["t₁"], p["t₂"], p["aF"], sss.ω_vec)
            ft_F = ft_pert_ramp(p["t₁"], p["t₂"], Fmax, a, sss.ω_vec, -get_c(x0)*Δx )

            get_ft_F = LinearInterpolation(sss.f_vec, abs.(ft_F))
            fres_sample = range(sss.f1_res, stop=sss.f2_res, length=sss.n_int)
            freq_product = get_ft_F.( fres_sample ) .* get_ampresp.( fres_sample )
            A_f_res_integrated = sum( freq_product ) / sss.n_int
            append!(resfreq_scat, A_f_res_integrated)
        end
    end
    # println( sum(resfreq_scat), sum(x_scat) )
    return [Fmax_scat, a_scat, resfreq_scat, x_scat]
end

function plot_scatter(x, sss, grid_axs, grid_fig, Δx)
    # println(sum(x[3]), sum(x[4]))
    scatter!(grid_axs[1], x[2], x[1], color = x[4], colormap = sss.cb_maps[1], colorrange = sss.cb_limits[1])
    scatter!(grid_axs[2], x[2], x[1], color = x[3], colormap = sss.cb_maps[2], colorrange = sss.cb_limits[2])
    save_fig("plots_v0.2/tmp/animations/", "Deltax=$Δx", "both", grid_fig)
end