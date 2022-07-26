#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using DifferentialEquations,
    LinearAlgebra, CairoMakie, Colors, Interpolations, Statistics, NLsolve;

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("ews_plo.jl"))
include(srcdir("plots_plo.jl"))
include(srcdir("video_helper_plo.jl"))
include(srcdir("mechanics_plo.jl"))
include(srcdir("pwlinear_oscillator.jl"))

#################################################
#################################################
#################################################

# Define plotting options.
plot_characteristics_bool = false
plot_bifportrait_bool = false
plot_superposition_bool = false
plot_tip_grid_bool = true

# Decide whether animation should be created or not.
plot_types = ["single", "Δx", "D", "σ"]             # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[3]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("tmp/"), plotsdir()]
prefix = prefixes[1]

p = load_highfreq_parameters()

#################################################
#################################################
#################################################

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
    p["fixed_tend"] = true
else
    p["noisy"] = false
    p["fixed_dt"] = false
    p["fixed_tend"] = false
end

p["Δx"] = 0.6                           # Standard value if Δx not changed over looping.-

# p["fixed_dt"] = true
# p["dt"] = 1e-1                      # accelerate computation

if plot_type == "single"
    title_func(x) = ""
elseif plot_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = [0, 0.3, 0.6, 0.9]

    # In case we want to generate an animation:
    gen_anim = false
    Δx_vec_anim = range(0, stop = 1.0, step = 0.1)    # Vector of sampled initial conditions.

elseif plot_type == "D"
    title_func(x) = "D = $x"
    # D_vec = [1e-1, 2e-1, 5e-1, 1]    # Vector of sampled dampings.
    D_vec = [0.0, 1e-4, 1e-3, 1e-2]    # Vector of sampled dampings.
elseif plot_type == "σ"
    title_func(x) = "σ = $x N"
    σ_vec = [0.25, 0.5]
end
p["F_type"] = "ramp"

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
# a_llim, a_ulim = -2, 3                  # Range of sampled slopes.
a_llim, a_ulim = -1, 3                  # Range of sampled slopes.
Δx = p["Δx"]                            # Standard value if Δx not changed over looping.-
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

sol_dict = Dict()
if plot_tip_grid_bool

    # Initialise axis and dictionary for tipping grid.
    grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
    grid_axs = init_grid_axs(grid_fig)
    grid_dict = Dict()

    if plot_type == "single"
        plot_scatter(
            get_scatter(Δx, plot_type, sss, solve_plo_F),
            sss,
            grid_axs,
            grid_fig,
            Δx,
            prefix_anim,
            "fixed_cb",
        )
    else
        if plot_type == "Δx"
            for Δx in Δx_vec
                grid_dict[string(Δx)] =
                    get_scatter(Δx, plot_type, sss, solve_plo_F)
            end
            plot_grid4(grid_dict, Δx_vec, prefix)

            if gen_anim
                record(
                    grid_fig,
                    string(prefix_anim, "slices.mp4"),
                    Δx_vec;
                    framerate = framerate,
                ) do Δx
                    grid_dict[string(Δx)] =
                        get_scatter(Δx, plot_type, sss, solve_plo_F)
                    plot_scatter(
                        grid_dict[string(Δx)],
                        sss,
                        grid_axs,
                        grid_fig,
                        Δx,
                        prefix_anim,
                    )
                end
            end
        elseif plot_type == "D"
            for D in D_vec
                grid_dict[string(D)] =
                    get_scatter(D, plot_type, sss, solve_plo_F)
            end
            plot_grid4(grid_dict, D_vec, prefix)

        elseif plot_type == "σ"
            for σ in σ_vec
                grid_dict[string(σ)] =
                    get_scatter(σ, plot_type, sss, solve_plo_F)
            end
            plot_grid2(grid_dict, σ_vec, prefix)
        end
    end
end