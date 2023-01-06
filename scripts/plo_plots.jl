using DrWatson
@quickactivate "nlo-rtipping"

using CairoMakie, JLD2, Colors
include(srcdir("utils.jl"))

files = [ "grid4_Δx.jld2", "grid4_D.jld2", "EM_grid_tbuffer10_dt0.01_nm100.jld2"]
gdicts = [ JLD2.load( datadir(files[i]), "grid_dict" ) for i in eachindex(files) ]

Δ_vals = [0.0, 0.6, 0.9]'
D_vals = [1e-2, 0.5, 1.0]'
σ_vals = [.3, 1., 2.]'
vals = vcat(Δ_vals, D_vals, σ_vals)

for i in eachindex(gdicts)
    for key in keys(gdicts[i])
        nkey = string( round( parse( Float64, key ); digits = 3 ) )
        gdicts[i][nkey] = gdicts[i][key]
    end
end

tks = ( 10. .^ (-2:3), [L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$", L"$10^{2}$", L"$10^{3}$"])
colorranges = [(1., 3.), (1., 3.), (0, 1)]
sep = [[2.25], [2.25], [0.49]]
cmaps = [:rainbow1, :rainbow1, :rainbow]

fig = Figure(resolution = (1500, 1500), font = srcdir("cmunrm.ttf"), fontsize = 28)
nrows, ncols = 3, 3

for i in 1:nrows, j in 1:ncols
    gdict = gdicts[i]
    local node = vals[i, j]
    local titles = [
        L"$\Delta x_1 = %$(string(node)) \, \mathrm{m}$",
        L"$D = %$(string(node))",
        L"$\sigma = %$(string(node)) \, \mathrm{N}$" ]

    x = gdict[string(node)]
    ax = Axis(
        fig[i, j][1, 1],
        title = titles[i],
        xlabel = (i == nrows ? L"$a \: \mathrm{(N \, s^{-1})}$" : " "),
        ylabel = (j == 1 ? L"$F_{\max}$ (N)" : " "),
        xscale = log10,
        xaxisposition = :bottom,
        yaxisposition = (j == 1 ? :left : :right),
        xticks = tks,
        xticklabelsvisible = true,
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
        colormap = cmaps[i],
        colorrange = colorranges[i],
        markersize = 5,
    )

    contour!(
        ax,
        x[2],
        x[1],
        x[3],
        color = :black,
        levels = sep[i],
        linewidth = 2,
    )

end

Colorbar(
    fig[1:2, ncols + 1],
    colormap = cmaps[1],
    colorrange = colorranges[1],
    label = L"$x_{1}(t = t_{e})$ (m)",
    height = Relative( .5 ),
    lowclip = cgrad(cmaps[1])[1],
    highclip = cgrad(cmaps[1])[end],
)

Colorbar(
    fig[3, ncols + 1],
    colormap = cmaps[3],
    colorrange = colorranges[3],
    label = L"$\hat{P}_\mathrm{tip}$ ",
    height = Relative( .8 ),
)

save_fig(plotsdir("plo/"), "plo_grid", "both", fig)