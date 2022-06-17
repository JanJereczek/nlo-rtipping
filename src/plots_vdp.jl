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
    fig = Figure( resolution = (1500, 750), font = srcdir("cmunrm.ttf"), fontsize=28 )
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

function plot_grid4(scatter_dict, node_vec, prefix)
    fig = Figure( resolution = (1500, 1500), font = srcdir("cmunrm.ttf"), fontsize=28 )
    pws = [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"]
    # title = L"$\Delta x_{1} =$ %$(string(Δx)) m",
    for i in 1:2
        for j in 1:2
            l = (i-1)*2 + j
            node = node_vec[l]
            x = scatter_dict[string(node)]
            ax = Axis(
                fig[i,j][1,1], 
                title = L"$\,$ %$(string(node))", 
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
            hm = scatter!(ax, x[2], x[1], color = x[3], colormap = :rainbow1, colorrange = (0, 3) )
            # ct = contour!(ax, x[2], x[1], x[3], color = :black, levels = [2.25, 2.26], linewidth = 5)
        end
    end
    Colorbar(fig[:, 3], colormap = :rainbow1, colorrange = (0, 3), height = Relative(1/2), highclip = :red, label = L"$|| x(t = t_{e}) ||$")
    save_fig(prefix, "grid4", "both", fig)
end


function plot_superposition(Fvec, avec, Δx, p)
    nrows = 2
    ncols = 3

    fig = Figure(resolution = (1500, 800), font = srcdir("cmunrm.ttf"), fontsize=28 )
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
            local sol = solve_plo(x₀, tspan, p)

            local solF_ = solve_plo_F([0., 0.], tspan, p)
            local solF = solF_(t)

            p["aF"] = 0
            local soly0_ = solve_plo(x₀, tspan, p)
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