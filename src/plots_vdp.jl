include("utils.jl")
include("pwlinear_oscillator.jl")

function init_grid_axs(fig)
    ax = Axis(
        fig[1, 1][1, 1],
        ylabel = L"$F_{\max}$ [N]",
        xlabel = L"$a$ [N/s]",
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

function plot_grid2(scatter_dict, node_vec, prefix)
    fig = Figure(resolution = (1500, 750), font = srcdir("cmunrm.ttf"), fontsize = 28)
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    for i = 1:2
        node = node_vec[i]
        x = scatter_dict[string(node)]
        ax = Axis(
            fig[1, i][1, 1],
            title = L"$\sigma =$ %$(string(node)) N",
            xlabel = L"$a$ [N/s]",
            ylabel = L"$F_{\max}$ [N]",
            xscale = log10,
            xticks = (10.0 .^ (-2:3), pws),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
        )
        hm = scatter!(
            ax,
            x[2],
            x[1],
            color = x[3],
            colormap = :rainbow1,
            colorrange = (1.25, 3.25),
        )
    end
    Colorbar(
        fig[:, 3],
        colormap = :rainbow1,
        colorrange = (1.25, 3.25),
        label = L"$x_{1}(t = t_{e})$ [m]",
    )
    save_fig(prefix, "grid2", "both", fig)
end

function plot_grid_vdp(scatter_dict, node_vec, prefix)
    fig = Figure(resolution = (1500, 500), font = srcdir("cmunrm.ttf"), fontsize = 28)
    pws = [L"$10^{-3}$", L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$"]
    for i in eachindex(node_vec)
        node = node_vec[i]
        x = scatter_dict[string(node)]
        titles = [
            L"$\Delta x = $ %$(string(node))",
            L"$\mu = $ %$(string(node))",
        ]
        ax = Axis(
            fig[1, i][1, 1],
            title = titles[selector],
            xlabel = L"$a$ ",
            ylabel = (i == 1 ? L"$F_{\max}$ " : " "),
            xscale = log10,
            xaxisposition = :bottom,
            yaxisposition = :left ,
            xticks = (10.0 .^ (-3:1), pws),
            yminorticks = IntervalsBetween(5),
            yminorgridvisible = true,
            yticklabelsvisible = (i == 1 ? true : false),
        )
        hm = scatter!(
            ax,
            x[2],
            x[1],
            color = x[3],
            colormap = :rainbow1,
            highclip = cgrad(:rainbow1)[end],
            lowclip = cgrad(:rainbow1)[1],
            colorrange = (p["F_crit"] + F_llim, 3),
            markersize = 5,
        )
        contour!(
            ax,
            x[2],
            x[1],
            x[3],
            color = :black,
            levels = [3],
            linewidth = 2,
        )
    end
    Colorbar(
        fig[:, 4],
        colormap = :rainbow1,
        colorrange = (p["F_crit"] + F_llim, 3),
        height = Relative( .8 ),
        highclip = cgrad(:rainbow1)[end],
        lowclip = cgrad(:rainbow1)[1],
        label = L"$|| x(t = t_{e}) ||$",
    )
    save_fig(prefix, "_grid", "both", fig)
end

function plot_superposition(
    Fvec,
    avec,
    Δx,
    p,
    solve_nlo,
    solve_nlo_F;
    Fapprox = 45,
    ia = 30,
)
    nrows = 2
    ncols = 3

    fig = Figure(resolution = (1500, 800), font = srcdir("cmunrm.ttf"), fontsize = 28)
    Fmax = Fvec[argmin(abs.(Fvec .- Fapprox))]
    x₀ = get_x₀(p, Δx)
    tend = 100
    tspan = (0, tend)
    t = range(tspan[1], stop = tspan[2], step = 0.1)

    for i = 1:nrows
        for j = 1:ncols
            k = (i - 1) * ncols + j
            p["aF"] = avec[ia+k]
            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + p["Fmax"] / p["aF"]
            Fplot = get_Fvec(t)
            local sol_ = solve_nlo_F(x₀, tspan, p)
            local sol = sol_(t)

            local solF_ = solve_nlo_F([0.0, 0.0], tspan, p)
            local solF = solF_(t)

            p["aF"] = 0
            local soly0_ = solve_nlo(x₀, tspan, p)
            local soly0 = soly0_(t)

            ax = Axis(
                fig[i, j],
                title = string("Experiment ", k),
                xlabel = (i == nrows ? L"$t$ [s]" : ""),
                ylabel = (j == 1 ? L"$x_{1}$ [m]" : ""),
                xminorticks = IntervalsBetween(5),
                xminorgridvisible = true,
                yminorticks = IntervalsBetween(5),
                yminorgridvisible = true,
            )
            # ax0 = Axis(
            #     fig[i, j],
            #     ylabel = (j == 3 ? L"$F$ [N]" : ""),
            #     ylabelcolor = :gray,
            #     yticklabelcolor = :gray,
            #     yaxisposition = :right,
            #     yticks = 0:0.5:2,
            #     xgridvisible = false,
            #     ygridvisible = false,
            # )
            # hidespines!(ax0)
            # hidexdecorations!(ax0)

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
            # l0 = lines!(ax0, t, Fplot, label = L"$F(t)$", color = :gray)
            l1 = lines!(ax, sol.t, y[:, 1], label = L"$x_{1}(t)$")
            l2 = lines!(ax, solF.t, yF[:, 1], label = L"$x_{1}^{F}(t)$")
            l3 = lines!(ax, soly0.t, yy0[:, 1], label = L"$x_{1}^{g}(t)$")
            l4 = lines!(ax, t, ycheck[:, 1], label = L"$x_{1}^{F}(t) + x_{1}^{g}(t)$")
            l5 = hlines!(ax, [-2, 2], color = :red, label = L"$x_{T} = 2$ m")
            lins = [l1, l2, l3, l4, l5]
            ylims!(ax, (-2.5, 2.5))

            if (i == nrows) & (j == ncols)
                Legend(
                    fig[:, 4],
                    [l1, l2, l3, l4, l5],
                    [
                        L"$x_{1}(t)$",
                        L"$x_{1}^{F}(t)$",
                        L"$x_{1}^{g}(t)$",
                        L"$x_{1}^{F}(t) + x_{1}^{g}(t)$",
                        L"$x_{T} = 2$ m",
                    ],
                )
            end
        end
    end
    save_fig(prefix, "superposition", "both", fig)
end

function plot_stream_vdp(Fconst_vec, x1bnds, x2bnds)
    fig = Figure(resolution = (1200, 1200), font = srcdir("cmunrm.ttf"), fontsize = 20)
    nrows, ncols = 2, 2
    axs = []
    for i in 1:nrows, j in 1:ncols
        k = j + (i-1)*ncols
        local ax = Axis(
            fig[i, j],
            xlabel = L"$x_{1}$ [m]",
            ylabel = L"$x_{2}$ [m/s]",
            title = string("F=", Fconst_vec[k]),
        )
        streamF(u) = forced_vdP_stream( u, Fconst_vec[k] )
        sp = streamplot!(
            ax,
            streamF,
            x1bnds[k][1] .. x1bnds[k][2],
            x2bnds[k][1] .. x2bnds[k][2],
            gridsize = (40, 40),
        )
    end
    save_fig(plotsdir("vdp/"), "stream", "both", fig)
end