# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, Statistics, DrWatson, NLsolve;

# Include self-written scripts.
include(srcdir("forcing.jl"))
include(srcdir("mechanics.jl"))
include(srcdir("nonlinear_oscillator.jl"))
include(srcdir("parameters.jl"))
include(srcdir("plots.jl"))
include(srcdir("utils.jl"))
include(srcdir("video_helper.jl"))
include(srcdir("ews.jl"))

# Parameter dictionnary of the nlo.
p = load_parameters()

# Define plotting options.
plot_characteristics_bool = true
plot_bifportrait_bool = true
plot_tip_grid_bool = true
plot_superposition_bool = true

# Decide whether animation should be created or not.
plot_types = ["none", "Δx", "D", "σ"]      # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[2]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("tmp/"), plotsdir()]
prefix = prefixes[1]

if plot_type == "none"
    prefix_anim = string(prefix, "Δx")
else
    prefix_anim = string(prefix, "animations/", plot_type)
end

# Control run time.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -2, 3                  # Range of sampled slopes.

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

Δx = 0.6
p["F_type"] = "ramp"
if (plot_type == "none") | (plot_type == "grid4")
    title_func(x) = ""
elseif plot_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = [0, 0.3, 0.6, 0.9]

    # In case we want to generate an animation:
    gen_anim = false
    Δx_vec_anim = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.
    
elseif plot_type == "D"
    title_func(x) = "D = $x"
    D_vec = range(0.1, stop = 2, step = 0.1)    # Vector of sampled dampings.
    d_vec = D2d(D_vec)
elseif plot_type == "σ"
    title_func(x) = "σ = $x N"
    σ_vec = [0.25, 0.5]
end

# If D > 1 the equation has a singularity.
if p["D"] < 1
    ω_vec = 10 .^ (range(-4, stop = 4, length = 8001))
    res_threshold = 1.1
    get_ω_ampresp, get_f_ampresp, amp_ω_resp, ω_res, f_res =
        get_resonance_characteristics(p, ω_vec, res_threshold)
end

# Plot the characteristics
if plot_characteristics_bool
    # Get characteristics of spring n°2.
    x₁_vec = range(0.0, stop = 2.0, length = 1000)
    c₂ = get_nl_stiffness.(x₁_vec)
    plot_characteristics(x₁_vec, c₂, ω_vec, amp_ω_resp, prefix, "ω")
end

# Get system bifurcation.
if plot_bifportrait_bool
    Fbif = range(0, stop = 200, length = 2001)
    x_bassin_bound = get_bassin_boundary(p, 0.0)
    plot_bifurcation_stream(Fbif, prefix, x_bassin_bound)
end

# Sampled values of the tipping grid.
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
n_int = 100
node = Observable(0.0)
title_node = lift(title_func, node)
sss = SlicedScatterStructs(avec, Fvec, ω_vec, ω_res, n_int)

# Plot the superposed solutions.
if plot_superposition_bool
    plot_superposition(Fvec, avec, Δx, p)
end

if plot_tip_grid_bool

    # Initialise axis and dictionary for tipping grid.
    grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
    grid_axs = init_grid_axs(grid_fig, title_node)
    grid_dict = Dict()

    if plot_type == "none"
        plot_scatter(get_scatter(Δx, plot_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim, "fixed_cb")
    else
        if plot_type == "Δx"
            for Δx in Δx_vec
                grid_dict[string(Δx)] = get_scatter(Δx, plot_type, sss)
            end
            plot_grid4(grid_dict, Δx_vec, prefix)

            if gen_anim
                record(grid_fig, string(prefix_anim, "slices.mp4"), Δx_vec; framerate = framerate) do Δx
                    node[] = Δx
                    plot_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim)
                end
            end
        elseif plot_type == "D"
            for D in D_vec
                grid_dict[string(D)] = get_scatter(D, plot_type, sss)
            end
            plot_grid2(grid_dict, D_vec, prefix)
        elseif plot_type == "σ"
            for σ in σ_vec
                grid_dict[string(σ)] = get_scatter(σ, plot_type, sss)
            end
            plot_grid2(grid_dict, σ_vec, prefix)
        end
    end
end