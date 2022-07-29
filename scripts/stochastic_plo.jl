#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using CairoMakie, Colors, ProgressMeter;
using DifferentialEquations;
using LinearAlgebra, Interpolations, Statistics, NLsolve;
using .Threads, BenchmarkTools, JLD2;

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
compute_grid_bool = false
plot_tip_grid_bool = true

plot_type = "σ"
prefix = plotsdir("plo/stochastic/")

p = load_parameters()
p["noisy"] = true
p["fixed_dt"] = true
p["dt"] = 1e-2
p["fixed_tend"] = true
p["t_buffer"] = 10
p["Δx"] = 0.6
p["F_type"] = "ramp"
σ_vec = [1, 2, 5, 10]

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

grid_dict = Dict()
nm = 100
grid_file = string("stochastic_grid_tbuffer", p["t_buffer"], "_dt", p["dt"], "_nm", nm, ".jld2")

if compute_grid_bool
    for σ in σ_vec
        grid_dict[string(σ)] = threaded_stochastic_scatter(σ, sss, solve_plo, nm)
    end
    jldsave(datadir(grid_file); grid_dict)
end

if plot_tip_grid_bool
    grid_dict = JLD2.load( datadir(grid_file), "grid_dict" )
    plot_stochastic_grid(grid_dict, σ_vec, prefix)
end