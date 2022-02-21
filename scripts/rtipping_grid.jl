# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, ProgressMeter, Statistics, DrWatson;

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
plot_tip_grid_bool = false
plot_response_bool = false

# Decide whether animation should be created or not.
anim_types = ["none", "Δx", "D", "σ"]   # Choose between slicing (or not) over IC, damping degree or noise variance.
anim_type = anim_types[1]
framerate = 2                          # [frame/second]

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
    title_func(x) = ""
    Δx = 0.6
elseif anim_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.
elseif anim_type == "D"
    title_func(x) = "D = $x"
    D_vec = range(0.1, stop = 2, step = 0.1)    # Vector of sampled dampings.
    d_vec = D2d(D_vec)
elseif anim_type == "σ"
    title_func(x) = "σ = $x N"
    σ_vec = range(0.0, stop = 1.0, step = 0.05)
end


ω_vec = 10 .^ (range(-4, stop = 4, length = 8001))
res_threshold = 1.1
get_ω_ampresp, get_f_ampresp, amp_ω_resp, ω_res, f_res =
    get_resonance_characteristics(p, ω_vec, res_threshold)

# Plot the characteristics
if plot_characteristics_bool
    # Get characteristics of spring n°2.
    x₁_vec = range(0.0, stop = 2.0, length = 1000)
    c₂ = get_nl_stiffness.(x₁_vec)
    plot_characteristics(x₁_vec, c₂, ω_vec, amp_ω_resp, prefix, "ω")
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

# Get the tranform of the seperating control trajectory.
region = 1
p["Fmax"] = 0.95*p["F_crit"]
Δx_lim = - ( p["Fmax"] ) / get_c(p["xeq1"])
tf_lim = Dict()
p["x₀_lim"] = get_x₀(p, Δx_lim)[1]
p["Ystat"] = get_Y_stat(p["x₀_lim"], ω_vec)

# tf_lim["U"] = ft_step( p["Fmax"] + p["m"]*p["g"] , ω_vec)
# tf_lim["G"], tf_lim["Y₀"] = nlo_transfer(p, tf_lim["x₀"], ω_vec, region)
# tf_lim["Y"] = get_Y(tf_lim["G"], tf_lim["Y₀"], tf_lim["U"])
# fig_lim = plot_bode(p, ω_vec, string(prefix, "lim_"), tf_lim)

# Sampled values of the tipping grid.
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)

n_int = 100
cb_maps = [:rainbow, :vik]
# cb_maps = [:rainbow, :thermal]
# cb_limits = [(1.4, 3.1), (415, 680)]
cb_limits = [(1.4, 3.1), (-0.2, 0.2)]

node = Observable(0.0)
title_node = lift(title_func, node)

# Initialise figure of tipping grid.
grid_fig = Figure(resolution = (1600, 800), fontsize = 18)
grid_axs = init_grid_axs(grid_fig, title_node)
sss = SlicedScatterStructs(avec, Fvec, ω_vec, ω_res, n_int, cb_limits, cb_maps)

function get_bode_a(a, p, Δx, Fvec, ω_vec)
    tf = Dict()
    tf["x₀"] = get_x₀(p, Δx)
    p_temp = p
    j = 28
    j = 32
    p_temp["Fmax"], p_temp["aF"], p_temp["t₂"] = Fvec[j], a, Fvec[j] / a

    tf["U"] = ft_stepramp(p_temp["t₁"], p_temp["t₂"], p_temp["m"]*p_temp["g"], p_temp["aF"], ω_vec)
    tf["G"], tf["Y₀"] = nlo_transfer(p_temp, tf["x₀"], ω_vec, region)
    tf["Y"] = get_Y(tf["G"], tf["Y₀"], tf["U"])
    tf["Ystat"] = p["Ystat"]
    println("a=", round(a; digits=2), ":  ", sum(abs.(tf["Y"]) .> abs.(tf_lim["Ystat"])))
    plot_bode(p, ω_vec, string(prefix, "a=", a), tf)
end

if anim_type == "none"
    plot_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim, "fixed_cb")
    # plot_scatter(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim)
    # for i in 30:40
    #     get_bode_a(avec[i], p, Δx, Fvec, ω_vec)
    # end
else
    Colorbar(grid_fig[1, 1][1, 2], colormap = cb_maps[1], limits = cb_limits[1], label = "x₁(t_end) [m]")
    Colorbar(grid_fig[1, 2][1, 2], colormap = cb_maps[2], limits = cb_limits[2], label = "|Y(ω_res)|")
    if anim_type == "Δx"
        record(grid_fig, string(prefix_anim, "slices.mp4"), Δx_vec; framerate = framerate) do Δx
            node[] = Δx
            plot_bode(p, Δx, ω_vec, prefix)
            plot_scatter_fixedcb(get_scatter(Δx, anim_type, sss), sss, grid_axs, grid_fig, Δx, prefix_anim)
        end
    elseif anim_type == "D"
        record(grid_fig, string(prefix_anim, "slices.mp4"), D_vec; framerate = framerate) do D
            node[] = D
            plot_scatter_fixedcb(get_scatter(D, anim_type, sss), sss, grid_axs, grid_fig, D, prefix_anim)
        end
    elseif anim_type == "σ"
        record(grid_fig, string(prefix_anim, "slices.mp4"), σ_vec; framerate = framerate) do σ
            node[] = σ
            plot_scatter_fixedcb(get_scatter(σ, anim_type, sss), sss, grid_axs, grid_fig, σ, prefix_anim)
        end
    end
end