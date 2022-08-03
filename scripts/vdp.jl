#################################################
############# Initialization ####################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using DifferentialEquations, CairoMakie, Colors
using Interpolations, Statistics, NLsolve, ProgressMeter

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("plots_vdp.jl"))
include(srcdir("integrate_tipgrid.jl"))
include(srcdir("vanderpol_oscillator.jl"))

# Define plotting options.
plot_stream = false
plot_superposition_bool = false
compute_grid_bool = true
plot_tip_grid_bool = true
plot_types = ["Δx", "μ"]
selector = 1
plot_type = plot_types[selector]

# Define saving options.
prefixes = [plotsdir("vdp/"), plotsdir()]
prefix = string( prefixes[1], plot_type )

# Set default values if parameter not varied.
p = load_parameters()
p["Δx"] = 1.0
p["μ"] = 0.05
Δx = p["Δx"]
μ = p["μ"]

# Set integration options.
p["noisy"] = false
p["fixed_dt"] = false
p["fixed_tend"] = true
p["F_type"] = "ramp"

# Define sampled values of initial perturbation and damping.
Δx_vec = [0, 0.8, 1.6]
μ_vec = [1e-2, 1e-1, 1e0]
node_vecs = [Δx_vec, μ_vec]
node_vec = node_vecs[selector]

#################################################
#################################################
#################################################

# Sampled values of the tipping grid.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -2, 1                  # Range of sampled slopes.
F_llim, F_ulim = -0.5, .2               # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
sss = SlicedScatterStructs(avec, Fvec)

if compute_grid_bool
    # Dictionaries to store solutions for tipping grid.
    grid_dict = Dict()
    sol_dict = Dict()
    for node in node_vec
        get_scatter!(node, plot_type, sss, solve_forced_vdP, grid_dict, sol_dict)
    end
    grid_file = string("vdp_", plot_type, ".jld2")
    jldsave(datadir(grid_file); grid_dict)
end

# if plot_tip_grid_bool
#     grid_dict = JLD2.load( datadir(grid_file), "grid_dict" )
#     plot_grid_vdp(grid_dict, node_vec, prefix)
# end

#################################################
############### Optional plots ##################
#################################################

if plot_stream
    Fconst_vec = [0, 1, 1.2, 1.5]
    x1bnds = [(-3, 3), (0.5, 1.5), (0.5, 1.5), (1, 2)]
    x2bnds = [(-3, 3), (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]
    plot_stream_vdp(Fconst_vec, x1bnds, x2bnds)
end

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