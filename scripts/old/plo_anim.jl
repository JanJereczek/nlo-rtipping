if gen_anim
    record(
        grid_fig,
        string(prefix_anim, "slices.mp4"),
        Δx_vec;
        framerate = framerate,
    ) do Δx
        grid_dict[string(Δx)] =
            get_scatter(Δx, plot_type, sss, solve_plo_F)
        plot_scatter(
            grid_dict[string(Δx)],
            sss,
            grid_axs,
            grid_fig,
            Δx,
            prefix_anim,
        )
    end
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