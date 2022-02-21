include("utils.jl")
include("nonlinear_oscillator.jl")

function plot_equil_control(p, prefix)
    ctrl_x₀ = [p["xₜ"] - 1e-1, 0.0]
    sol = solve_nlo(ctrl_x₀, [0, 100], p)
    u = reduce(vcat, transpose.(sol.u))

    fig = Figure(resolution = (400, 800), fontsize = 14)
    ax1 = Axis(fig[1, 1], ylabel = "F [N]")
    ax2 = Axis(fig[2, 1], xlabel = "Time [s]", ylabel = "x [m]")

    lines!(ax1, sol.t, get_F.(sol.t))
    lines!(ax2, sol.t, u[:, 1])
    ylims!(ax1, (100, 150))
    ylims!(ax2, (0, 2))
    save_fig(prefix, "control_run", "both", fig)
end

function plot_characteristics(x₁_vec, c₂, f_vec, amp_resp, prefix, mode)
    fig = Figure(resolution = (1000, 500), fontsize = 14)

    if mode == "ω"
        freqlabel = "f [rad/s]"
    else
        freqlabel = "f [Hz]"
    end

    ax1 = Axis(
        fig[1, 1],
        title = "Characteristic curve of spring n°2",
        xlabel = "x [m]",
        ylabel = "c₂(x) [kg/s²]",
    )

    ax2 = Axis(
        fig[1, 2],
        title = "Amplitude response before tipping",
        xlabel = freqlabel,
        ylabel = "Amplification [dB]",
        xscale = log10,
        yscale = log10,
        xminorticks = IntervalsBetween(9),
        xminorticksvisible = true,
        xminorgridvisible = true,
    )

    lines!(ax1, x₁_vec, c₂)
    lines!(ax2, f_vec, amp_resp)
    save_fig(prefix, "characteristics", "both", fig)
end

function plot_bifurcation(Fbif, prefix)
    fig = Figure(resolution = (400, 400), fontsize = 14)
    ax = Axis(
        fig[1, 1],
        title = "Bifurcation diagramm of the system",
        xlabel = "F + mg [N]",
        ylabel = "x [m]",
    )
    lines!(ax, Fbif, branch_xeq1.(Fbif))
    lines!(ax, Fbif, branch_xeq2.(Fbif))
    save_fig(prefix, "bifurcation", "both", fig)
end

function plot_phaseportrait(X1, X2, prefix)
    X = cat(vec(X1), vec(X2), dims = 2)
    fig = Figure(resolution = (800, 400), fontsize = 14)
    ax1 = Axis(fig[1, 1], title = "Vector field", xlabel = "x₁ [m]", ylabel = "x₂ [m/s]")
    ax2 = Axis(fig[1, 2], title = "Phase portrait", xlabel = "x₁ [m]", ylabel = "x₂ [m/s]")

    dX = mapslices(nl_osc_free, X; dims = 2)
    mag = vec(sqrt.(dX[:, 1] .^ 2 .+ dX[:, 2] .^ 2))

    arrows!(
        ax1,
        X1,
        X2,
        dX[:, 1],
        dX[:, 2],
        arrowsize = 3,
        lengthscale = 0.05,
        arrowcolor = mag,
        linecolor = mag,
    )
    streamplot!(
        ax2,
        nl_osc_free_stream,
        0 .. 3,
        -5 .. 5,
        colormap = :rainbow1,
        gridsize = (32, 32),
        arrow_size = 0,
    )
    save_fig(prefix, "phaseportrait", "both", fig)
end

function init_grid_axs(fig, title_node)
    ax1 = Axis(
        fig[1, 1][1, 1],
        title = title_node,
        ylabel = "Fₘ [N]",
        xlabel = "a [N/s]",
        xscale = log10,
    )
    ax2 = Axis(
        fig[1, 2][1, 1],
        title = "Spectral amplitude at resonance frequency",
        ylabel = "Fₘ [N]",
        xlabel = "a [N/s]",
        xscale = log10,
    )
    axs = [ax1, ax2]

    return axs
end

function init_large_grid_axs(fig, title_node, nrows, ncols)
    axs = Dict()
    for i in 1:nrows
        for j in 1:ncols
            k = (i-1)*ncols + j
            if k == nrows*ncols     # init time-integrated tipping grid.
                ax = Axis(
                    fig[i,j][1, 1],
                    title = title_node,
                    xlabel = "a [N/s]",
                    ylabel = "Fₘ [N]",
                    xscale = log10,
                )
            else                    # init the measure plots.
                ax = Axis(fig[i,j][1,1], title=string("γ", k), xlabel = "a [N/s]", ylabel = "Fₘ [N]", xscale = log10)
            end
            axs[string(k)] = ax
        end
    end
    return axs
end

function init_response_axs(fig, f1, f2)
    ax1 = Axis(fig[1, 1], xlabel = "Time [s]", ylabel = "F(t) [N]")
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Frequency [Hz]",
        ylabel = "Amplitude [1]",
        xscale = log10,
        xminorticks = IntervalsBetween(9),
        xminorticksvisible = true,
        xminorgridvisible = true,
    )
    ax3 = Axis(fig[2, 1], xlabel = "Time [s]", ylabel = "x₁(t) [m]")
    ax4 = Axis(fig[2, 2], xlabel = "x₁ [m]", ylabel = "x₂ [m/s]")
    axs = [ax1, ax2, ax3, ax4]
    lines!(axs[2], [f1, f1], [-10, 2000], color = :red)
    lines!(axs[2], [f2, f2], [-10, 2000], color = :red)
    return axs
end

# function plot_spect_grid(fig, ax, fmax_scat, a_scat, resfreq_scat)

# end

function show_response(x₀, Fmax, a, t, p, f_vec, ft, axs)
    F = get_Fvec(t)
    sol = solve_nlo(x₀, (t[1], t[end]), p)
    line_label = string("Fmax=", round(Fmax, digits = 2), ",  a=", round(a, digits = 2))
    lines!(axs[1], t, F)
    lines!(axs[2], f_vec, ft)
    lines!(axs[3], sol.t, sol[1, :])
    lines!(axs[4], sol[1, :], sol[2, :], label = line_label)
    # ylims!(axs[2])
end

function plot_bode(p::Dict, ω_vec::Vector, prefix::String, tf::Dict)

    G, Y₀, U, Y, Ystat = tf["G"], tf["Y₀"], tf["U"], tf["Y"], tf["Ystat"]
    fig = Figure(resolution = (1000, 800))
    ax1 = Axis(
        fig[1, 1],
        xlabel = "ω [rad/s]",
        ylabel = "Magnitude [dB]",
        xscale = log10,
        yscale = log10,
        xticks =  10 .^ range(-2., 1., step=1),
        xminorticks = IntervalsBetween(9),
        xminorticksvisible = true,
        xminorgridvisible = true,
    )
    ax2 = Axis(
        fig[2, 1],
        xlabel = "ω [rad/s]",
        ylabel = "Phase [rad]",
        xscale = log10,
        xticks =  10 .^ range(-2., 1., step=1),
        xminorticks = IntervalsBetween(9),
        xminorticksvisible = true,
        xminorgridvisible = true,
    )
    xlims!(ax1, (5e-2, 1e1))
    xlims!(ax2, (5e-2, 1e1))

    tf_vec = [G, Y₀, U, Y, Ystat]
    label_vec = ["G", "Y₀", "U", "Y", "Ystat"]

    for i = 1:length(tf_vec)
        lines!(ax1, ω_vec, abs.(tf_vec[i]), label = label_vec[i])
        lines!(ax2, ω_vec, angle.(tf_vec[i]) )
    end

    vlines!(ax1, [p["ω₀1"]], color = :red)
    vlines!(ax2, [p["ω₀1"]], color = :red)
    ylims!(ax1, (1e-3, 1e3))
    ylims!(ax2, (-4, .5))

    axislegend(ax1, merge = true, nbanks = 2, position = :rt)
    save_fig(prefix, "Bodeplot_Δx=$Δx", "both", fig)
    return fig
end
