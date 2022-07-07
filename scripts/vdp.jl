#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, Statistics, NLsolve;

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("plots_vdp.jl"))
include(srcdir("video_helper.jl"))
include(srcdir("vanderpol_oscillator.jl"))

#################################################
#################################################
#################################################

# Define plotting options.
plot_stream = false
plot_superposition_bool = true
plot_tip_grid_bool = false

# Decide whether animation should be created or not.
plot_types = ["single", "Δx", "D", "σ"]             # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[2]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("vdp/"), plotsdir()]
prefix = prefixes[1]

#################################################
#################################################
#################################################

p = load_parameters()
p["μ"] = 0.02

if plot_stream
    fig = Figure(resolution = (800, 800), font = srcdir("cmunrm.ttf"), fontsize = 20)
    ax = Axis(fig[1, 1][1, 1], xlabel = L"$x_{1}$ [m]", ylabel = L"$x_{2}$ [m/s]")
    sp = streamplot!(ax, vdP_stream, -3 .. 3, -3 .. 3, gridsize = (32, 32))
    save_fig(plotsdir("vdp/"), "stream", "both", fig)
end

#################################################
#################################################
#################################################

if plot_type == "single"
    prefix_anim = string(prefix, "Δx")
else
    prefix_anim = string(prefix, "animations/", plot_type)
end

# Choose applied noise.
if plot_type == "σ"
    p["noisy"] = true
    prefix = string(prefix, "noisy_")
    p["fixed_dt"] = true
    p["dt"] = 1e-1                      # [s] SDE requires fixed time step (or more advanced concepts).
    p["fixed_tend"] = true              # should be fixed in order to have same energy input regardless of time scale.
else
    p["noisy"] = false
    p["fixed_dt"] = false
    p["fixed_tend"] = true
end

if plot_type == "single"
    title_func(x) = ""
elseif plot_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = [0, 0.5, 1.0, 1.5]
    # Δx_vec = [1.5, 1.8, 2., 2.5]

    # In case we want to generate an animation:
    gen_anim = false
    Δx_vec_anim = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.

elseif plot_type == "D"
    title_func(x) = "D = $x"
    D_vec = [1e-2, 5e-2, 1e-1, 1.5]    # Vector of sampled dampings.
elseif plot_type == "σ"
    title_func(x) = "σ = $x N"
    σ_vec = [0.25, 0.5]
end
p["F_type"] = "ramp"

#################################################
#################################################
#################################################

# Sampled values of the tipping grid.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -3, 1                  # Range of sampled slopes.
p["Δx"] = 0.6                           # Standard value if Δx not changed over looping.-
Δx = p["Δx"]                           # Standard value if Δx not changed over looping.-
F_llim, F_ulim = -1.5, -0.8                 # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
sss = SlicedScatterStructs(avec, Fvec)

# Plot the superposed solutions.
if plot_superposition_bool
    plot_superposition(
        Fvec,
        avec,
        1.5,
        p,
        solve_vdP,
        solve_forced_vdP;
        Fapprox = 0.7,
        ia = 25,
    )
end

#################################################
#################################################
#################################################

sol_dict = Dict()
if plot_tip_grid_bool

    # Initialise axis and dictionary for tipping grid.
    grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
    grid_axs = init_grid_axs(grid_fig)
    grid_dict = Dict()

    if plot_type == "single"
        plot_scatter(
            get_scatter(Δx, plot_type, sss, solve_forced_vdP),
            sss,
            grid_axs,
            grid_fig,
            Δx,
            prefix_anim,
            "fixed_cb",
        )
    else
        if plot_type == "Δx"
            for Δx in Δx_vec
                grid_dict[string(Δx)], sol_dict[string(Δx)] =
                    get_scatter(Δx, plot_type, sss, solve_forced_vdP)
            end
            plot_grid4(grid_dict, Δx_vec, prefix)

            if gen_anim
                record(
                    grid_fig,
                    string(prefix_anim, "slices.mp4"),
                    Δx_vec;
                    framerate = framerate,
                ) do Δx
                    grid_dict[string(Δx)], sol_dict[string(Δx)] =
                        get_scatter(Δx, plot_type, sss, solve_forced_vdP)
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
        elseif plot_type == "D"
            for D in D_vec
                grid_dict[string(D)], sol_dict[string(Δx)] =
                    get_scatter(D, plot_type, sss, solve_forced_vdP)
            end
            plot_grid4(grid_dict, D_vec, prefix)

        elseif plot_type == "σ"
            for σ in σ_vec
                grid_dict[string(σ)], sol_dict[string(Δx)] =
                    get_scatter(σ, plot_type, sss, solve_forced_vdP)
            end
            plot_grid2(grid_dict, σ_vec, prefix)
        end
    end
end