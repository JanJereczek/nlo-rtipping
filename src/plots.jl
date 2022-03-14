include("utils.jl")
include("nonlinear_oscillator.jl")

function plot_equil_control(p, prefix)
    ctrl_x₀ = [p["xₜ"] - 1e-1, 0.0]
    sol = solve_nlo(ctrl_x₀, [0, 100], p)
    u = reduce(vcat, transpose.(sol.u))

    fig = Figure(resolution = (400, 800), fontsize = 14)
    ax1 = Axis(fig[1, 1], ylabel = L"$F(t)$ [N]")
    ax2 = Axis(fig[2, 1], xlabel = L"$t$ [s]", ylabel = L"$x$ [m]")

    lines!(ax1, sol.t, get_F.(sol.t))
    lines!(ax2, sol.t, u[:, 1])
    ylims!(ax1, (100, 150))
    ylims!(ax2, (0, 2))
    save_fig(prefix, "control_run", "both", fig)
end

function plot_characteristics(x₁_vec, c₂, f_vec, amp_resp, prefix, mode)
    fig = Figure(resolution = (1000, 500), fontsize = 14)

    if mode == "ω"
        freqlabel = L"$\omega$ [rad/s]"
    else
        freqlabel = L"$f$ [Hz]"
    end

    ax1 = Axis(
        fig[1, 1],
        title = "Characteristic curve of spring n°2",
        xlabel = L"$x$ [m]",
        ylabel = L"$c_{2}(x)$ [kg/s²]",
    )

    ax2 = Axis(
        fig[1, 2],
        xlabel = freqlabel,
        ylabel = L"$|\cdot|$ [dB]",
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

function plot_bifurcation_stream(Fbif, prefix, x_bassin_bound)
    fig = Figure(resolution = (1200, 500), font = "/home/jan/pCloudDrive/My Documents/Fonts/cmu/cmunrm.ttf", fontsize = 20)
    ax1 = Axis(
        fig[1, 1],
        xlabel = L"$\tilde{F} + mg$ [N]",
        ylabel = L"$\tilde{x}$ [m]",
    )
    lines!(ax1, Fbif, branch_xeq1.(Fbif), label = L"$\tilde{x}_{-}$ [m]")
    lines!(ax1, Fbif, branch_xeq2.(Fbif), label = L"$\tilde{x}_{+}$ [m]")
    axislegend(ax1, merge = true, nbanks = 2, position = :lt)
    
    ax2 = Axis(fig[1, 2][1, 1], xlabel = L"$x_{1}$ [m]", ylabel = L"$x_{2}$ [m/s]")

    streamplot!(
        ax2,
        nl_osc_free_stream,
        -1 .. 3,
        -8 .. 5,
        gridsize = (32, 32),
    )
    lines!(ax2, x_bassin_bound[1,:], x_bassin_bound[2,:], linewidth = 3, color = :red)
    lines!(ax2, [x_bassin_bound[1, 1], x_bassin_bound[1, end]], [x_bassin_bound[2, 1], x_bassin_bound[2, end]], linewidth = 3, color = :red )
    save_fig(prefix, "bif_phaseportrait", "both", fig)
end

function init_grid_axs_1(fig)
    ax = Axis(
        fig[1, 1][1, 1],
        ylabel = L"$F_{\max}$ [N]",
        xlabel = L"$a$ [N/s]",
        xscale = log10,
    )
    return ax
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

function show_response(x₀, Fmax, a, t, p, f_vec, ft, axs)
    F = get_Fvec(t)
    sol = solve_nlo(x₀, (t[1], t[end]), p)
    line_label = string("Fmax=", round(Fmax, digits = 2), ",  a=", round(a, digits = 2))
    lines!(axs[1], t, F)
    lines!(axs[2], f_vec, ft)
    lines!(axs[3], sol.t, sol[1, :])
    lines!(axs[4], sol[1, :], sol[2, :], label = line_label)
end

function plot_grid2(scatter_dict, node_vec, prefix)
    fig = Figure( resolution = (1500, 750), font = "/home/jan/pCloudDrive/My Documents/Fonts/cmu/cmunrm.ttf", fontsize=28 )
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    for i in 1:2
        node = node_vec[i]
        x = scatter_dict[string(node)]
        ax = Axis(
            fig[1,i][1,1], 
            title = L"$\sigma =$ %$(string(node)) N", 
            xlabel = L"$a$ [N/s]",
            ylabel = L"$F_{\max}$ [N]",
            xscale = log10,
            xticks = (10. .^ (-2:3), pws),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
        )
        hm = scatter!(ax, x[2], x[1], color = x[3], colormap = :rainbow1, colorrange = (1.25, 3.25) )
    end
    Colorbar(fig[:, 3], colormap = :rainbow1, colorrange = (1.25, 3.25), label = L"$x_{1}(t = t_{e})$ [m]")
    save_fig(prefix, "grid2", "both", fig)
end

function plot_grid4(scatter_dict, Δx_vec, prefix)
    fig = Figure( resolution = (1500, 1500), font = "/home/jan/pCloudDrive/My Documents/Fonts/cmu/cmunrm.ttf", fontsize=28 )
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    for i in 1:2
        for j in 1:2
            l = (i-1)*2 + j
            Δx = Δx_vec[l]
            x = scatter_dict[string(Δx)]
            ax = Axis(
                fig[i,j][1,1], 
                title = L"$\Delta x =$ %$(string(Δx)) m", 
                xlabel = L"$a$ [N/s]",
                # xlabel = (i == 1 ? L"$\hat{t}$ [s]" : L"$a$ [N/s]"),
                ylabel = L"$F_{\max}$ [N]",
                xscale = log10,
                xaxisposition = (i == 1 ? :top : :bottom),
                yaxisposition = (j == 1 ? :left : :right),
                xticks = (10. .^ (-2:3), pws),
                yminorticks = IntervalsBetween(5),
                yminorgridvisible = true,
            )
            hm = scatter!(ax, x[2], x[1], color = x[3], colormap = :rainbow1, colorrange = (1.25, 3.25) )
            ct = contour!(ax, x[2], x[1], x[3], color = :black, levels = [2.25, 2.26], linewidth = 5)
        end
    end
    Colorbar(fig[:, 3], colormap = :rainbow1, colorrange = (1.25, 3.25), label = L"$x_{1}(t = t_{e})$ [m]")
    save_fig(prefix, "grid4", "both", fig)
end

function plot_superposition(Fvec, avec, Δx, p)
    nrows = 2
    ncols = 3

    fig = Figure(resolution = (1500, 800), font = "/home/jan/pCloudDrive/My Documents/Fonts/cmu/cmunrm.ttf", fontsize=28 )
    ia = 31
    Fmax = Fvec[ isapprox.(Fvec, 45; atol=1e-1) ][1]
    x₀ = get_x₀(p, Δx)
    tend = 10
    tspan = (0, tend)
    t = range(tspan[1], stop=tspan[2], step=0.1)
    
    for i in 1:nrows
        for j in 1:ncols
            k = (i-1)*ncols + j
            p["aF"] = avec[ ia+k ]
            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + p["Fmax"] / p["aF"]
            Fplot = get_Fvec(t)
            local sol = solve_nlo(x₀, tspan, p)

            local solF_ = solve_nlo_F([0., 0.], tspan, p)
            local solF = solF_(t)

            p["aF"] = 0
            local soly0_ = solve_nlo(x₀, tspan, p)
            local soly0 = soly0_(t)

            ax = Axis(
                    fig[i, j],
                    title = string("Experiment ", k),
                    xlabel = (i == nrows ? L"$t$ [s]" : ""),
                    ylabel = (j == 1 ? L"$x_{1}$ [m]" : ""), 
                    xminorticks=IntervalsBetween(5), 
                    xminorgridvisible=true,
                    yminorticks=IntervalsBetween(5), 
                    yminorgridvisible=true,
                )
            ax0 = Axis(fig[i, j], ylabel = (j == 3 ? L"$F$ [N]" : ""), ylabelcolor = :gray, yticklabelcolor = :gray, yaxisposition = :right, yticks = 0:50:150, xgridvisible = false, ygridvisible = false)
            hidespines!(ax0)
            hidexdecorations!(ax0)

            if j < ncols
                hideydecorations!(ax0)
            end
            if j > 1
                hideydecorations!(ax, grid = false, minorgrid = false)
            end
            if i < nrows
                hidexdecorations!(ax, grid = false, minorgrid = false)
            end

            y = reduce(vcat, transpose.(sol.u))
            yF = reduce(vcat, transpose.(solF.u))
            yy0 = reduce(vcat, transpose.(soly0.u))
            ycheck = yF .+ yy0
            l0 = lines!(ax0, t, Fplot, label=L"$F(t)$", color = :gray)
            l1 = lines!(ax, sol.t, y[:, 1], label=L"$x_{1}(t)$")
            l2 = lines!(ax, solF.t, yF[:, 1], label=L"$x_{1}^{F}(t)$")
            l3 = lines!(ax, soly0.t, yy0[:, 1], label=L"$x_{1}^{g}(t)$")
            l4 = lines!(ax, t, ycheck[:, 1], label=L"$x_{1}^{F}(t) + x_{1}^{g}(t)$")
            l5 = hlines!(ax, [p["xₜ"]], color = :red, label = L"$x_{T} = 1.5$ m")
            lins = [l0, l1, l2, l3, l4, l5]
            ylims!(ax, (0, 1.6))
            ylims!(ax0, (0, 160))
            if ( i == nrows ) & ( j == ncols )
                Legend(fig[:, 4], [l0, l1, l2, l3, l4, l5],
                       [L"$F(t)$", L"$x_{1}(t)$", L"$x_{1}^{F}(t)$", L"$x_{1}^{g}(t)$", L"$x_{1}^{F}(t) + x_{1}^{g}(t)$", L"$x_{T} = 1.5$ m"])
            end
        end
    end
    save_fig(prefix, "superposition", "both", fig)
end