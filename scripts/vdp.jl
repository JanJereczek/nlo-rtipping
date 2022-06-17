# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, Statistics, DrWatson, NLsolve;

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("mechanics.jl"))
include(srcdir("plots_vdp.jl"))
include(srcdir("video_helper_vdp.jl"))
include(srcdir("vanderpol_oscillator.jl"))

# Define plotting options.
plot_superposition_bool = false
plot_tip_grid_bool = true

# Decide whether animation should be created or not.
plot_types = ["single", "Δx", "D", "σ"]             # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[2]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("vdp/"), plotsdir()]
prefix = prefixes[1]

p = load_parameters()

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
else
    p["noisy"] = false
    p["fixed_dt"] = false
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

# Sampled values of the tipping grid.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -2, 3                  # Range of sampled slopes.
Δx = 0.3                                # Standard value if Δx not changed over looping.-
F_llim, F_ulim = -1, .5                 # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
sss = SlicedScatterStructs(avec, Fvec)

# Plot the superposed solutions.
if plot_superposition_bool
    plot_superposition(Fvec, avec, Δx, p)
end

sol_dict = Dict()
if plot_tip_grid_bool

    # Initialise axis and dictionary for tipping grid.
    grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
    grid_axs = init_grid_axs(grid_fig)
    grid_dict = Dict()

    if plot_type == "single"
        plot_scatter(get_scatter(Δx, plot_type, sss, solve_forced_vdP), sss, grid_axs, grid_fig, Δx, prefix_anim, "fixed_cb")
    else
        if plot_type == "Δx"
            for Δx in Δx_vec
                grid_dict[string(Δx)], sol_dict[string(Δx)] = get_scatter(Δx, plot_type, sss, solve_forced_vdP)
            end
            plot_grid4(grid_dict, Δx_vec, prefix)

            if gen_anim
                record(grid_fig, string(prefix_anim, "slices.mp4"), Δx_vec; framerate = framerate) do Δx
                    grid_dict[string(Δx)], sol_dict[string(Δx)] = get_scatter(Δx, plot_type, sss, solve_forced_vdP)
                    plot_scatter(grid_dict[string(Δx)], sss, grid_axs, grid_fig, Δx, prefix_anim)
                end
            end
        elseif plot_type == "D"
            for D in D_vec
                grid_dict[string(D)], sol_dict[string(Δx)] = get_scatter(D, plot_type, sss, solve_forced_vdP)
            end
            plot_grid4(grid_dict, D_vec, prefix)

        elseif plot_type == "σ"
            for σ in σ_vec
                grid_dict[string(σ)], sol_dict[string(Δx)] = get_scatter(σ, plot_type, sss, solve_forced_vdP)
            end
            plot_grid2(grid_dict, σ_vec, prefix)
        end
    end
end