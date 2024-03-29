include("utils.jl")
include("pwlinear_oscillator.jl")

function plot_equil_control(p, prefix)
    ctrl_x₀ = [p["xₜ"] - 1e-1, 0.0]
    sol = solve_plo(ctrl_x₀, [0, 100], p)
    u = reduce(vcat, transpose.(sol.u))

    fig = Figure(resolution = (400, 800), fontsize = 14)
    ax1 = Axis(fig[1, 1], ylabel = L"$F(t)$ (N)")
    ax2 = Axis(fig[2, 1], xlabel = L"$t$ (s)", ylabel = L"$x$ (m)")

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
        xlabel = L"$x$ (m)",
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
    fig = Figure(resolution = (1200, 500), font = srcdir("cmunrm.ttf"), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = L"$\tilde{F} + mg$ (N)", ylabel = L"$\tilde{x}$ (m)", title = L"(a) $\,$")
    lines!(ax1, 0.0:0.1:150, branch_xeq1.(0.0:0.1:150), label = L"$\tilde{x}_{-}$ (m)")
    lines!(ax1, Fbif, branch_xeq2.(Fbif), label = L"$\tilde{x}_{+}$ (m)", linestyle = :dash)
    axislegend(ax1, merge = true, nbanks = 2, position = :lt)

    ax2 = Axis(fig[1, 2], xlabel = L"$x_{1}$ (m)", ylabel = L"$x_{2}$ (m/s)", title = L"(b) $\,$")
    ax3 = Axis(fig[1, 3], xlabel = L"$x_{1}$ (m)", yticklabelsvisible = false, title = L"(c) $\,$")

    crange = (0, 10)
    cmap = cgrad([:grey20, :grey20]) #:jet
    sp = streamplot!(
        ax2,
        pwlin_oscillator_stream,
        -1 .. 3, -5 .. 5,
        gridsize = (32, 32),
        colormap = cmap,
        colorrange = crange,
    )
    k1_tmp = copy(p["k₁"])
    p["k₁"] = 0.0
    sp = streamplot!(
        ax3,
        pwlin_oscillator_stream,
        -1 .. 3, -5 .. 5,
        gridsize = (32, 32),
        colormap = cmap,
        colorrange = crange,
    )
    p["k₁"] = k1_tmp

    # Colorbar(
    #     fig[2,2:3],
    #     colormap = cmap,
    #     colorrange = crange,
    #     width = Relative(0.5),
    #     label = L"Normalised 1-norm of $\dot{x}$ (1)",
    #     vertical = false,
    # )

    # streamplot!(ax2, pwlin_oscillator_stream, -1 .. 3, -8 .. 5, gridsize = (32, 32))
    # lines!(ax2, x_bassin_bound[1, :], x_bassin_bound[2, :], linewidth = 3, color = :red)
    # lines!(
    #     ax2,
    #     [x_bassin_bound[1, 1], x_bassin_bound[1, end]],
    #     [x_bassin_bound[2, 1], x_bassin_bound[2, end]],
    #     linewidth = 3,
    #     color = :red,
    # )
    save_fig(prefix, "bif_phaseportrait", "both", fig)
end

function init_grid_axs(fig)
    ax = Axis(
        fig[1, 1][1, 1],
        ylabel = L"$F_{\max}$ (N)",
        xlabel = L"$a$ (N/s)",
        xscale = log10,
    )
    return ax
end

function show_response(x₀, Fmax, a, t, p, f_vec, ft, axs)
    F = get_Fvec(t)
    sol = solve_plo(x₀, (t[1], t[end]), p)
    line_label = string("Fmax=", round(Fmax, digits = 2), ",  a=", round(a, digits = 2))
    lines!(axs[1], t, F)
    lines!(axs[2], f_vec, ft)
    lines!(axs[3], sol.t, sol[1, :])
    lines!(axs[4], sol[1, :], sol[2, :], label = line_label)
end

function plot_stochastic_grid(scatter_dict, node_vec, prefix)
    fig = Figure(resolution = (1500, 700), font = srcdir("cmunrm.ttf"), fontsize = 28)
    pws = [L"$10^{0}$", L"$10^{1}$", L"$10^{2}$"]
    nrows, ncols = 1, 3
    for i in 1:nrows, j in 1:ncols
        k = (i - 1) * ncols + j
        node = node_vec[k+1]
        x = scatter_dict[string(node)]
        ax = Axis(
            fig[i, j][1, 1],
            title = L"$\sigma =$ %$(string(node)) N",
            xlabel = (i == nrows ? L"$a$ (N/s)" : " "),
            ylabel = (j == 1 ? L"$F_{\max}$ (N)" : " "),
            xscale = log10,
            # xaxisposition = (i == 1 ? :top : :bottom),
            yaxisposition = (j == 1 ? :left : :right),
            xticks = (10.0 .^ (0:2), pws),
            # xticklabelsvisible = (i == nrows ? true : false),
            yticklabelsvisible = (j == 1 ? true : false),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
            xminorgridvisible = true,
        )

        scatter_bool = true
        if scatter_bool
            hm = scatter!(
                ax,
                x[2],
                x[1],
                color = x[3],
                colormap = :rainbow1,
                colorrange = (0, 1),
                markersize = 5,
            )

            contour!(
                ax,
                x[2],
                x[1],
                x[3],
                color = :black,
                levels = [0.5],
                linewidth = 2,
            )
        else
            hm = contourf!(
                ax,
                x[2],
                x[1],
                x[3],
                colormap = :rainbow1,
                levels = -.1:.1:1.1,
            )
        end
    end
    
    Colorbar(
        fig[:, ncols + 1],
        colormap = :rainbow1,
        colorrange = (0, 1),
        label = L"$\hat{P}_\mathrm{tip}$",
        height = Relative( .8 ),
    )
    save_fig(prefix, "stochastic_grid", "both", fig)
end

function plot_grid4(scatter_dict, node_vec, prefix)
    fig = Figure(resolution = (1500, 1200), font = srcdir("cmunrm.ttf"), fontsize = 28)
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    nrows, ncols = 2, 2
    clims = (1.25, 3.25)
    for i in 1:nrows, j in 1:ncols
        k = (i - 1) * ncols + j
        node = node_vec[k]
        x = scatter_dict[string(node)]
        rnode = round(node; digits = 3 )
        title_dict = Dict( 
            "Δx" => L"$\Delta x_{1} = %$(string( rnode )) \, $m",
            "D"  => L"$D = %$(string( rnode )) $",
            "σ"  => L"$\sigma = %$(string( rnode )) \,$N",
        )
        ax = Axis(
            fig[i, j][1, 1],
            title = title_dict[plot_type],
            xlabel = (i == nrows ? L"$a$ (N/s)" : " "),
            ylabel = (j == 1 ? L"$F_{\max}$ (N)" : " "),
            xscale = log10,
            xaxisposition = (i == 1 ? :top : :bottom),
            yaxisposition = (j == 1 ? :left : :right),
            xticks = (10.0 .^ (-2:3), pws),
            xticklabelsvisible = (i == nrows ? true : false),
            yticklabelsvisible = (j == 1 ? true : false),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
            xminorgridvisible = true,
        )

        hm = scatter!(
            ax,
            x[2],
            x[1],
            color = x[3],
            colormap = :rainbow1,
            colorrange = clims,
            markersize = 5,
        )
        ct = contour!(
            ax,
            x[2],
            x[1],
            x[3],
            color = :black,
            levels = [2.25, 2.26],
            linewidth = 2,
        )
    end
    Colorbar(
        fig[:, ncols + 1],
        colormap = :rainbow1,
        colorrange = clims,
        label = L"$x_{1}(t = t_{e})$ (m)",
        height = Relative( .5 ),
        highclip = :red,
        lowclip = :purple,
    )

    save_fig(prefix, string( plot_type, "_grid4" ), "both", fig)
end


function plot_grid3(scatter_dict, node_vec, prefix)
    fig = Figure(resolution = (1500, 700), font = srcdir("cmunrm.ttf"), fontsize = 28)
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    nrows, ncols = 1, 3
    clims = (1.25, 3.25)
    for i in 1:nrows, j in 1:ncols
        k = (i - 1) * ncols + j
        node = node_vec[k]
        x = scatter_dict[string(node)]
        rnode = round(node; digits = 3 )
        title_dict = Dict( 
            "Δx" => L"$\Delta x_{1} = %$(string( rnode )) \, $m",
            "D"  => L"$D = %$(string( rnode )) $",
            "σ"  => L"$\sigma = %$(string( rnode )) \,$N",
        )
        ax = Axis(
            fig[i+1, j][1, 1],
            title = title_dict[plot_type],
            xlabel = (i == nrows ? L"$a$ (N/s)" : " "),
            ylabel = (j == 1 ? L"$F_{\max}$ (N)" : " "),
            xscale = log10,
            yaxisposition = (j == 1 ? :left : :right),
            xticks = (10.0 .^ (-2:3), pws),
            xticklabelsvisible = (i == nrows ? true : false),
            yticklabelsvisible = (j == 1 ? true : false),
            yticksvisible = (j == 2 ? false : true),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
            xminorgridvisible = true,
        )

        hm = scatter!(
            ax,
            x[2],
            x[1],
            color = x[3],
            colormap = :rainbow1,
            colorrange = clims,
            markersize = 5,
        )
        ct = contour!(
            ax,
            x[2],
            x[1],
            x[3],
            color = :black,
            levels = [2.25, 2.26],
            linewidth = 2,
        )
    end
    Colorbar(
        fig[1, :],
        colormap = :rainbow1,
        colorrange = clims,
        label = L"$x_{1}(t = t_{e})$ (m)",
        highclip = :red,
        lowclip = :purple,
        vertical = false,
        width = Relative( .5 ),
    )

    save_fig(prefix, string( plot_type, "_grid4" ), "both", fig)
end


function plot_superposition(Fvec, avec, Δx, p, solve_nlo, solve_nlo_F; Fapprox = 45)
    nrows = 2
    ncols = 3

    fig = Figure(resolution = (1500, 900), font = srcdir("cmunrm.ttf"), fontsize = 28)
    ia = 31
    Fmax = Fvec[isapprox.(Fvec, Fapprox; atol = 1e-1)][1]
    x₀ = get_x₀(p, Δx)
    tend = 10
    tspan = (0, tend)
    t = range(tspan[1], stop = tspan[2], step = 0.01)

    for i = 1:nrows
        for j = 1:ncols
            k = (i - 1) * ncols + j
            p["aF"] = avec[ia+k]
            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + p["Fmax"] / p["aF"]
            Fplot = get_Fvec(t)
            local sol = solve_plo(x₀, tspan, p)
            local sol = sol(t)

            local solF_ = solve_plo_F([0.0, 0.0], tspan, p)
            local solF = solF_(t)

            p["aF"] = 0
            local soly0_ = solve_plo(x₀, tspan, p)
            local soly0 = soly0_(t)

            ax = Axis(
                fig[i, j],
                title = string("Experiment ", k),
                xlabel = (i == nrows ? L"$t$ (s)" : ""),
                ylabel = (j == 1 ? L"$x_{1}$ (m)" : ""),
                xminorticks = IntervalsBetween(5),
                xminorgridvisible = true,
                yminorticks = IntervalsBetween(5),
                yminorgridvisible = true,
            )
            ax0 = Axis(
                fig[i, j],
                ylabel = (j == 3 ? L"$F$ (N)" : ""),
                ylabelcolor = :gray,
                yticklabelcolor = :gray,
                yaxisposition = :right,
                yticks = 0:50:150,
                xgridvisible = false,
                ygridvisible = false,
            )
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
            l0 = lines!(ax0, t, Fplot, label = L"$F(t)$", color = :gray)
            l1 = lines!(ax, sol.t, y[:, 1], label = L"$x_{1}(t)$", color = :gray10)
            l2 = lines!(ax, solF.t, yF[:, 1], label = L"$x_{1}^{F}(t)$")
            l3 = lines!(ax, soly0.t, yy0[:, 1], label = L"$x_{1}^{g}(t)$", color = :green)
            l4 = lines!(ax, t, ycheck[:, 1], linestyle = :dash, label = L"$x_{1}^{F}(t) + x_{1}^{g}(t)$", color = :orange, linewidth = 2)
            l5 = hlines!(ax, [p["xₜ"]], color = :red, label = L"$x_{T} = 1.5$ m")
            lins = [l0, l1, l2, l3, l4, l5]
            ylims!(ax, (0, 1.6))
            ylims!(ax0, (0, 160))
            if (i == nrows) & (j == ncols)
                Legend(
                    fig[3, :],
                    [l0, l1, l2, l3, l4, l5],
                    [
                        L"$F(t)$",
                        L"$x_{1}(t)$",
                        L"$x_{1}^{F}(t)$",
                        L"$x_{1}^{g}(t)$",
                        L"$x_{1}^{F}(t) + x_{1}^{g}(t)$",
                        L"$x_{T} = 1.5$ m",
                    ],
                    orientation = :horizontal,
                    width = Relative(.8),
                    colgap = 80,
                )
            end
        end
    end
    save_fig(prefix, "superposition", "both", fig)
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