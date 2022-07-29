#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using DifferentialEquations, LinearAlgebra, NLsolve
using CairoMakie, Colors, Interpolations, Statistics, JLD2

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("plots_plo.jl"))
include(srcdir("mechanics_plo.jl"))
include(srcdir("pwlinear_oscillator.jl"))
include(srcdir("integrate_tipgrid_plo.jl"))

#################################################
#################################################
#################################################

# Define plotting options.
plot_characteristics_bool = false
plot_bifportrait_bool = false
plot_superposition_bool = false
plot_tip_grid_bool = true

# Decide whether animation should be created or not.
selector = 1                        # Selector for case difference.
plot_types = ["Δx", "D"]            # Slicing over IC or damping ratio...
plot_type = plot_types[selector]    # Choose which of both.
framerate = 2                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("tmp/"), plotsdir()]
prefix = prefixes[1]

p = ( selector == 1 ? load_parameters() : load_highfreq_parameters() )
p["noisy"] = false
p["fixed_dt"] = false
p["fixed_tend"] = false
p["F_type"] = "ramp"
p["Δx"] = 0.0                           # Standard value if Δx not changed over looping.
Δx = p["Δx"]
println(p["D"])

#################################################
#################################################
#################################################

prefix_anim = string(prefix, "animations/", plot_type)

if plot_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = [0, 0.3, 0.6, 0.9]

    # In case we want to generate an animation:
    gen_anim = false
    Δx_vec_anim = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.

elseif plot_type == "D"
    title_func(x) = "D = $x"
    # D_vec = [1e-1, 2e-1, 5e-1, 1]    # Vector of sampled dampings.
    D_vec = [1f-2, 1f-1, 5f-1, 1f0]    # Vector of sampled dampings.
end

# If D > 1 the equation has a singularity.
if p["D"] < 1
    ω_vec = 10 .^ (range(-4, stop = 4, length = 8001))
    res_threshold = 1.1
    get_ω_ampresp, get_f_ampresp, amp_ω_resp, ω_res, f_res =
        get_resonance_characteristics(p, ω_vec, res_threshold)
end

#################################################
#################################################
#################################################

# Plot the characteristics
if plot_characteristics_bool
    # Get characteristics of spring n°2.
    x₁_vec = range(0.0, stop = 2.0, length = 1000)
    c₂ = get_nl_stiffness.(x₁_vec)
    plot_characteristics(x₁_vec, c₂, ω_vec, amp_ω_resp, prefix, "ω")
end

#################################################
#################################################
#################################################

# Get system bifurcation.
if plot_bifportrait_bool
    Fbif = range(0, stop = 200, length = 2001)
    x_bassin_bound = get_bassin_boundary(p, 0.0)
    plot_bifurcation_stream(Fbif, prefix, x_bassin_bound)
end

#################################################
#################################################
#################################################

# Sampled values of the tipping grid.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
# a_llim, a_ulim = -2, 3                 # Range of sampled slopes.
a_llim, a_ulim = -2, 3                  # Range of sampled slopes.
F_llim, F_ulim = -15, 2                 # linscale with ref = F_crit
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
node_vecs = [Δx_vec, D_vec]
node_vec = node_vecs[selector]

if plot_tip_grid_bool
    # Dictionaries to store solutions for tipping grid.
    grid_dict = Dict()
    sol_dict = Dict()
    for node in node_vec
        get_scatter!(node, plot_type, sss, solve_plo, grid_dict, sol_dict)
    end
    plot_grid4(grid_dict, node_vec, prefix)
end

grid_file = string("grid4_", plot_type)
jldsave(datadir(grid_file); grid_dict)