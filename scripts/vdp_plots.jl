using DrWatson
@quickactivate "nlo-rtipping"

using CairoMakie, JLD2, Colors
include(srcdir("utils.jl"))

files = [ "vdp_Δx.jld2", "vdp_μ.jld2"]
gdicts = [ JLD2.load( datadir(files[i]), "grid_dict" ) for i in eachindex(files) ]

Δ_vals = [0, 0.9, 1.8]'
μ_vals = [5e-1, 1e0, 5e0]'
vals = vcat(Δ_vals, μ_vals)

for i in eachindex(gdicts)
    for key in keys(gdicts[i])
        nkey = string( round( parse( Float64, key ); digits = 3 ) )
        gdicts[i][nkey] = gdicts[i][key]
    end
end

tks = ( 10. .^ (-3:1), [L"$10^{-3}$", L"$10^{-2}$", L"$10^{-1}$", L"$10^{0}$", L"$10^{1}$"])
colorranges = [(.5, 3.), (.5, 3.)]
sep = [[6], [6]]
cmaps = [:rainbow1, :rainbow1]

fig = Figure(resolution = (1500, 1000), font = srcdir("cmunrm.ttf"), fontsize = 28)
nrows, ncols = 2, 3

for i in 1:nrows, j in 1:ncols
    gdict = gdicts[i]
    local node = vals[i, j]
    local titles = [
        L"$\Delta x = %$(string(node))$",
        L"$\mu = %$(string(node))"]

    x = gdict[string(node)]
    ax = Axis(
        fig[i, j][1, 1],
        title = titles[i],
        xlabel = (i == nrows ? L"$a$" : " "),
        ylabel = (j == 1 ? L"$F_{\max}$" : " "),
        xscale = log10,
        xaxisposition = :bottom,
        yaxisposition = (j == 1 ? :left : :right),
        xticks = tks,
        xticklabelsvisible = ( i == 2 ? true : false ),
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
    label = L"$|| x(t = t_{e}) ||$",
    height = Relative( .5 ),
    lowclip = cgrad(cmaps[1])[1],
    highclip = cgrad(cmaps[1])[end],
)

save_fig(plotsdir("vdp/"), "vdp_grid", "both", fig)