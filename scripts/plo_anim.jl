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