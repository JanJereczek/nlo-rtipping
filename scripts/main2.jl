# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, ProgressMeter, DrWatson;

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
plot_characteristics_bool = true
plot_bifurcation_bool = true
plot_stream_bool = true
plot_tip_grid_bool = false
plot_response_bool = false

# Decide whether animation should be created or not.
anim_types = ["none", "Δx", "D", "σ"]   # Choose between slicing (or not) over IC, damping degree or noise variance.
anim_type = anim_types[1]
framerate = 10                          # [frame/second]

# Define saving options.
prefixes = [plotsdir("tmp/"), plotsdir()]
prefix = prefixes[1]
prefix_anim = string(prefix, "animations/", anim_type)

# Control run time.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -2, 3                  # Range of sampled slopes.

# Choose applied noise.
noisy_forcing = false
if noisy_forcing
    p["F_type"] = "noisy"
    p["σ"] = 0.5
    prefix = string(prefix, "noisy_")
    p["dt"] = 1e-1                      # [s] SDE requires fixed time step (or more advanced concepts).
else
    p["F_type"] = "ramp"
end

if anim_type == "none"
    title_func(x) = "x₁(tₑ)"
elseif anim_type == "Δx"
    title_func(x) = "x₁(tₑ) for Δx = $x m"
    Δx_vec = range(0, stop = 1.0, step = 0.02)    # Vector of sampled initial conditions.
elseif anim_type == "D"
    title_func(x) = "x₁(tₑ) for D = $x"
    D_vec = range(0.1, stop = 2, step = 0.1)    # Vector of sampled dampings.
    d_vec = D2d(D_vec)
elseif anim_type == "σ"
    title_func(x) = "x₁(tₑ) for σ = $x N"
    σ_vec = range(0.0, stop = 1.0, step = 0.05)
end


ω_vec = 10 .^ (range(-3.0, stop = 2.0, length = 1001))
res_threshold = 1.1
get_ω_ampresp, get_f_ampresp, amp_resp, ω_res, f_res =
    get_resonance_characteristics(p, ω_vec, res_threshold)

# Plot the characteristics
if plot_characteristics_bool
    # Get characteristics of spring n°2.
    x₁_vec = range(0.0, stop = 2.0, length = 1000)
    c₂ = get_nl_stiffness.(x₁_vec)
    plot_characteristics(x₁_vec, c₂, ω_vec, amp_resp, prefix, "ω")
end

# Get system bifurcation.
if plot_bifurcation_bool
    Fbif = range(0, stop = 200, length = 201)
    plot_bifurcation(Fbif, prefix)
end

# Phase portrait and vector field.
if plot_stream_bool
    n₁, n₂ = 20, 20
    x1 = range(0, stop = 3, length = n₁)
    x2 = range(-5, stop = 5, length = n₂)
    X1, X2 = meshgrid(x1, x2)
    plot_phaseportrait(X1, X2, prefix)
end

# Sampled values of the tipping grid.
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)

n_int = 100
cb_maps = [:rainbow, :viridis]
cb_limits = [(1.5, 3.1), (50.0, 100.0)]

node = Observable(0.0)
title_node = lift(title_func, node)

# Initialise figure of tipping grid.
grid_fig = Figure(resolution = (1600, 800), fontsize = 18)
grid_axs = init_grid_axs(grid_fig, title_node)
sss = SlicedScatterStructs(avec, Fvec, ω_vec, ω_res, n_int, cb_limits, cb_maps)

Colorbar(grid_fig[1, 2], colormap = cb_maps[1], limits = cb_limits[1])
Colorbar(grid_fig[1, 4], colormap = cb_maps[2], limits = cb_limits[2])

# record(grid_fig, string(prefix, "ICslices.mp4"), Δx_vec; framerate = framerate) do Δx
#     init_pert[] = Δx
#     plot_scatter(get_scatter(Δx, sss), sss, grid_axs, grid_fig, Δx)
# end