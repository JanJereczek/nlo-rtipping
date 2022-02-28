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

# Parameter dictionnary of the nlo.
p = load_parameters()

# Define plotting options.
plot_characteristics_bool = false
plot_bifurcation_bool = false
plot_stream_bool = false
plot_tip_grid_bool = true
plot_superposition_bool = false

# Decide whether animation should be created or not.
anim_types = ["none", "grid4", "Δx", "D", "σ"]      # Choose between slicing (or not) over IC, damping degree or noise variance.
anim_type = anim_types[2]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("tmp/"), plotsdir()]
prefix = prefixes[1]

if anim_type == "none"
    prefix_anim = string(prefix, "Δx")
else
    prefix_anim = string(prefix, "animations/", anim_type)
end

# Control run time.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -2, 3                  # Range of sampled slopes.

# Choose applied noise.
if anim_type == "σ"
    p["F_type"] = "noisy"
    prefix = string(prefix, "noisy_")
    p["dt"] = 1e-1                      # [s] SDE requires fixed time step (or more advanced concepts).
else
    p["F_type"] = "ramp"
end

if (anim_type == "none") | (anim_type == "grid4")
    title_func(x) = ""
    Δx = 0.6
    Δx_vec = [0, 0.3, 0.6, 0.9]
elseif anim_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.
elseif anim_type == "D"
    title_func(x) = "D = $x"
    D_vec = range(0.1, stop = 2, step = 0.1)    # Vector of sampled dampings.
    d_vec = D2d(D_vec)
elseif anim_type == "σ"
    title_func(x) = "σ = $x N"
    σ_vec = [0.25, 0.5]
    Δx = 0.6
end

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
if plot_bifurcation_bool
    Fbif = range(0, stop = 200, length = 2001)
    plot_bifurcation_stream(Fbif, prefix)
end

# Phase portrait and vector field.
if plot_stream_bool
    n₁, n₂ = 20, 20
    x1 = range(0, stop = 3, length = n₁)
    x2 = range(-5, stop = 5, length = n₂)
    X1, X2 = meshgrid(x1, x2)
    plot_phaseportrait(X1, X2, prefix)
end

# Get the tranform of the seperating control trajectory.
region = 1
p["Fmax"] = 0.95*p["F_crit"]
Δx_lim = - ( p["Fmax"] ) / get_c(p["xeq1"])
tf_lim = Dict()
p["x₀_lim"] = get_x₀(p, Δx_lim)[1]
p["Ystat"] = get_Y_stat(p["x₀_lim"], ω_vec)

# Sampled values of the tipping grid.
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)

n_int = 100
cb_maps = [:rainbow, :rainbow, :rainbow, :rainbow, :rainbow, :rainbow]
cb_limits = ["none", (0, 1), (-0.1, 0.1), "none", "none", "none"]

node = Observable(0.0)
title_node = lift(title_func, node)
sss = SlicedScatterStructs(avec, Fvec, ω_vec, ω_res, n_int, cb_limits, cb_maps)

if plot_superposition_bool
    plot_superposition(Fvec, avec, Δx, p)
end

# Initialise figure of tipping grid.
# grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
# grid_axs = init_grid_axs(grid_fig, title_node)
# grid_axs = init_large_grid_axs(grid_fig, title_node, nrows, ncols)
grid_dict = Dict()
if plot_tip_grid_bool
    if anim_type == "none"
        # large_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim, nrows, ncols)
        plot_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim, "fixed_cb")
        # for i in 30:40
        #     get_bode_a(avec[i], p, Δx, Fvec, ω_vec)
        # end
    else
        # Colorbar(grid_fig[1, 1][1, 2], colormap = cb_maps[1], limits = cb_limits[1], label = "x₁(t_end) [m]")
        # Colorbar(grid_fig[1, 2][1, 2], colormap = cb_maps[2], limits = cb_limits[2], label = "|Y(ω_res)|")
        if anim_type == "Δx"
            record(grid_fig, string(prefix_anim, "slices.mp4"), Δx_vec; framerate = framerate) do Δx
                node[] = Δx
                # plot_bode(p, Δx, ω_vec, prefix)
                plot_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim)
            end
        elseif anim_type == "grid4"
            for Δx in Δx_vec
                grid_dict[string(Δx)] = get_scatter(Δx, anim_type, sss)
            end
            plot_grid4(grid_dict, Δx_vec, prefix)
        elseif anim_type == "D"
            for D in D_vec
                grid_dict[string(D)] = get_scatter(D, anim_type, sss)
            end
            plot_grid2(grid_dict, D_vec, prefix)
        elseif anim_type == "σ"
            for σ in σ_vec
                grid_dict[string(σ)] = get_scatter(σ, anim_type, sss)
            end
            plot_grid2(grid_dict, σ_vec, prefix)
        end
    end
end