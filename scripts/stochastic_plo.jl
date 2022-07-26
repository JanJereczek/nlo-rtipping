#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using CairoMakie, Colors, ProgressMeter;
using DifferentialEquations;
using LinearAlgebra, Interpolations, Statistics, NLsolve;
using .Threads, BenchmarkTools;

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("plots_plo.jl"))
include(srcdir("video_helper.jl"))
include(srcdir("mechanics_plo.jl"))
include(srcdir("pwlinear_oscillator.jl"))
include(srcdir("stochastic_ensemble.jl"))

#################################################
#################################################
#################################################

# Define plotting options.
plot_characteristics_bool = false
plot_bifportrait_bool = false
plot_superposition_bool = false
plot_tip_grid_bool = true

plot_type = "σ"
prefix = plotsdir("plo/stochastic/")

p = load_parameters()
p["noisy"] = true
p["fixed_dt"] = true
p["dt"] = 1e-1
p["fixed_tend"] = true
p["Δx"] = 0.6
p["F_type"] = "ramp"

σ_vec = [0.25, 0.5]
# σ_vec = [0., 0.]
title_func(x) = "σ = $x N"

#################################################
#################################################
#################################################

# Sampled values of the tipping grid.
nF, na = 20, 20                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = 0, 2                   # Range of sampled slopes.
Δx = p["Δx"]                            # Standard value if Δx not changed over looping.-
F_llim, F_ulim = -7, 2                  # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
sss = SlicedScatterStructs(avec, Fvec)

# Plot the superposed solutions.
if plot_superposition_bool
    plot_superposition(Fvec, avec, Δx, p, solve_plo, solve_plo_F; Fapprox = Fvec[end])
end

#################################################
#################################################
#################################################

sol_dict = Dict()
grid_dict = Dict()
grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
grid_axs = init_grid_axs(grid_fig)

for σ in σ_vec
    grid_dict[string(σ)] = threaded_stochastic_scatter(σ, sss, solve_plo, 100)
end

plot_grid2(grid_dict, σ_vec, prefix)