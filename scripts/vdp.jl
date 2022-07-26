#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using DifferentialEquations, CairoMakie, Colors, Interpolations, Statistics, NLsolve;

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("forcing.jl"))
include(srcdir("plots_vdp.jl"))
include(srcdir("video_helper.jl"))
include(srcdir("vanderpol_oscillator.jl"))

#################################################
#################################################
#################################################

# Define plotting options.
plot_stream = false
plot_superposition_bool = false
plot_tip_grid_bool = true

# Decide whether animation should be created or not.
plot_types = ["Δx", "μ"]             # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[1]

# Define saving options.
prefixes = [plotsdir("vdp/"), plotsdir()]
prefix = prefixes[1]

#################################################
#################################################
#################################################

p = load_parameters()
# p["μ"] = 0.02
p["μ"] = 0.2

if plot_stream
    fig = Figure(resolution = (1200, 1200), font = srcdir("cmunrm.ttf"), fontsize = 20)
    nrows, ncols = 2, 2
    axs = []
    Fconst_vec = [0, 1, 1.2, 1.5]
    x1bnds = [(-3, 3), (0.5, 1.5), (0.5, 1.5), (1, 2)]
    x2bnds = [(-3, 3), (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]
    for i in 1:nrows, j in 1:ncols
        k = j + (i-1)*ncols
        local ax = Axis(
            fig[i, j],
            xlabel = L"$x_{1}$ [m]",
            ylabel = L"$x_{2}$ [m/s]",
            title = string("F=", Fconst_vec[k]),
        )
        streamF(u) = forced_vdP_stream( u, Fconst_vec[k] )
        sp = streamplot!(
            ax,
            streamF,
            x1bnds[k][1] .. x1bnds[k][2],
            x2bnds[k][1] .. x2bnds[k][2],
            gridsize = (40, 40),
        )
        # append!(axs, ax)
    end
    save_fig(plotsdir("vdp/"), "stream", "both", fig)
end

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
    p["fixed_tend"] = true              # should be fixed in order to have same energy input regardless of time scale.
else
    p["noisy"] = false
    p["fixed_dt"] = false
    p["fixed_tend"] = true
end

if plot_type == "Δx"
    title_func(x) = "Δx = $x m"
    Δx_vec = [0, 0.5, 1.0, 1.5]
elseif plot_type == "D"
    title_func(x) = "D = $x"
    D_vec = [1e-2, 5e-2, 1e-1, 1.5]    # Vector of sampled dampings.
end
p["F_type"] = "ramp"

#################################################
#################################################
#################################################

# Sampled values of the tipping grid.
nF, na = 50, 50                         # Number of points sampled in the ramp-parameter space.
a_llim, a_ulim = -3, 1                  # Range of sampled slopes.
p["Δx"] = 0.6                           # Standard value if Δx not changed over looping.-
Δx = p["Δx"]                            # Standard value if Δx not changed over looping.-
F_llim, F_ulim = -0.5, .2               # linscale with ref = F_crit
Fvec = round.(p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits = 5)
avec = round.(10 .^ (range(a_llim, stop = a_ulim, length = na)); digits = 5)
sss = SlicedScatterStructs(avec, Fvec)

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

#################################################
#################################################
#################################################

sol_dict = Dict()
if plot_tip_grid_bool

    # Initialise axis and dictionary for tipping grid.
    grid_fig = Figure(resolution = (1500, 1000), fontsize = 18)
    grid_axs = init_grid_axs(grid_fig)
    grid_dict = Dict()

    if plot_type == "Δx"
        for Δx in Δx_vec
            grid_dict[string(Δx)], sol_dict[string(Δx)] =
                get_scatter(Δx, plot_type, sss, solve_forced_vdP)
        end
        plot_grid4(grid_dict, Δx_vec, prefix)
    elseif plot_type == "μ"
        for D in D_vec
            grid_dict[string(D)], sol_dict[string(Δx)] =
                get_scatter(D, plot_type, sss, solve_forced_vdP)
        end
        plot_grid4(grid_dict, D_vec, prefix)

    end
end