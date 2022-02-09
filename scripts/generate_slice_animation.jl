# Import packages.
using DifferentialEquations, CairoMakie, Colors, ProgressMeter, Interpolations;

# Include self-written scripts.
include("forcing.jl")
include("mechanics.jl")
include("nonlinear_oscillator.jl")
include("parameters.jl")
include("plots.jl")
include("utils.jl")
include("video_helper.jl")

# Define saving prefix.
prefix = "plots_v0.2/tmp/animations/d=70"

# Define options.
gen_Δx_video = true
gen_D_video = false

# Control video generation run time.
framerate = 10           # frame/second
nF, na = 50, 50
Δx_vec = range(0, stop=1, step=.02)
a_llim, a_ulim = -2, 3      # logscale with ref = 10^0

# Parameter dictionnary.
p = load_parameters()
p["t₁"] = 0.0         # [s] time at which the ramp begins.
p["x_eq"] = 1.0
p["d"] = 70
p["F_type"] = "ramp"

# Get frequency characterstics.
f_res, ω_res = round.(get_resonance_freq(p); digits = 3)

ω_vec = 10 .^ (range(-3.0, stop = 2.0, length = 1001))
f_vec = ω_vec ./ (2 * π)
amp_resp = get_bode_amp.(ω_vec)
get_ampresp = LinearInterpolation(f_vec, amp_resp)

fres_vec = f_vec[amp_resp .> 1.1]
f1_res, f2_res = round(fres_vec[1]; digits = 3), round(fres_vec[end]; digits = 3)
println("Resonance freq (Hz) - lower bound, peak and upper bound: ", f1_res, "  ", f_res, "  ", f2_res)


# Sampled values of the tipping grid.
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
Fvec = round.( p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nf); digits=5)
avec = round.( 10 .^ (range(a_llim, stop = a_ulim, length = na)); digits=5)

n_int = 100
cb_maps = [ :rainbow, :viridis ]
cb_limits = [(1.5, 3.1), (50.0, 100.0)]

# Initialise figure of tipping grid.
init_pert = Node(0.0)
title_func(x) = "x₁(tₑ) for Δx = $x m"
title_node = lift(title_func, init_pert)
grid_fig = Figure(resolution = (1600, 800), fontsize = 18)

ax1 = Axis(
    grid_fig[1, 1],
    title = title_node,
    ylabel = "Fₘ [N]",
    xlabel = "a [N/s]",
    xscale = log10,
)
ax2 = Axis(
    grid_fig[1, 3],
    title = "Integrated spectral product around resonance",
    ylabel = "Fₘ [N]",
    xlabel = "a [N/s]",
    xscale = log10,
)
grid_axs = [ax1, ax2]

sss = SlicedScatterStructs(avec, Fvec, ω_vec, f_vec, f1_res, f2_res, n_int, cb_limits, cb_maps)

Colorbar(grid_fig[1, 2], colormap = cb_maps[1], limits = cb_limits[1])
Colorbar(grid_fig[1, 4], colormap = cb_maps[2], limits = cb_limits[2])

record(grid_fig, string(prefix, "ICslices.mp4"), Δx_vec; framerate = framerate) do Δx
    init_pert[] = Δx
    plot_scatter(get_scatter(Δx, sss), sss, grid_axs, grid_fig, Δx)
end