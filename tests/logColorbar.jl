using CairoMakie 

x = 10.0.^(1:0.1:4)
y = 1.0:0.1:5.0
data = x .* ones(Float64, 1, length(y))
fig = Figure()
cmap = cgrad(:viridis, scale=:log10)
ax, hm = heatmap(fig[1, 1], x, y, data; colormap=cmap,
    axis=(;xscale=log10,
            xminorticksvisible=true,
            xminorticks=IntervalsBetween(9))
            )
cb = Colorbar(fig[1, 2], hm;
    minorticksvisible=true,
    minorticks=IntervalsBetween(9),
    scale=log10
)
fig