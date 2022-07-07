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

# Define plotting options.
plot_superposition_bool = false
plot_tip_grid_bool = true

# Decide whether animation should be created or not.
plot_types = ["single", "Δx", "D", "σ"]             # Choose between slicing (or not) over IC, damping degree or noise variance.
plot_type = plot_types[2]
framerate = 2                                       # [frame/second]

# Define saving options.
prefixes = [plotsdir("vdp/"), plotsdir()]
prefix = prefixes[1]

p = load_parameters()
p["fixed_dt"] = false
p["μ"] = 4

fig = Figure(resolution = (800, 800), font = srcdir("cmunrm.ttf"), fontsize = 20)
ax = Axis(fig[1, 1][1, 1],
    xlabel = L"$x_{1}$ [m]",
    ylabel = L"$x_{2}$ [m/s]",
)
sp = streamplot!(
    ax,
    vdP_stream,
    -3 .. 3,
    -3 .. 3,
    gridsize = (32, 32),
)

case = 5

if case == 1
    ϵ = 1e-2
    x0_vec = [[-2 + ϵ, 0], [2 - ϵ, 0], [-2 - ϵ, 0], [2 + ϵ, 0], [0, 2], [0, -2]]
    for x0 in x0_vec
        sol = solve_vdP(x0, (0, 100), p)
        xsol = [ sol.u[i][j] for i in 1:length(sol.u), j in 1:2]
        lines!(ax, xsol[:, 1], xsol[:, 2], label = string("x0 = ", x0), linewidth = 3)
    end

elseif case == 2
    F0_vec = 0:.1:1
    for F0 in F0_vec
        p["F₀"] = F0
        sol = solve_constantly_forced_vdP([0., 0.], (0, 100), p)
        xsol = [ sol.u[i][j] for i in 1:length(sol.u), j in 1:2]
        lines!(ax, xsol[:, 1], xsol[:, 2], label = string("F0 = ", F0), linewidth = 3)
    end

elseif case == 3
    F = 1.5
    a_vec = 10. .^ collect(-3:1:2)
    p["F_type"] = "ramp"
    p["noisy"] = false
    for a in a_vec
        p["aF"] = a
        # tend = F / minimum(a_vec) + 10
        tend = F / a + 10
        tspan = (0.0, tend)
        p["t₂"] = p["t₁"] + F / a
        sol = solve_forced_vdP([0., 0.], tspan, p)
        xsol = [ sol.u[i][j] for i in 1:length(sol.u), j in 1:2]
        lines!(ax, xsol[:, 1], xsol[:, 2], label = string("a = ", a), linewidth = 3)
    end
elseif case == 4
    fig = Figure(resolution = (1600, 300), font = srcdir("cmunrm.ttf"), fontsize = 20)
    F0_vec = 0:.2:1
    axs = [Axis( fig[1,j], xlabel = L"$x_{1}$ [m]", ylabel = L"$x_{2}$ [m/s]", ) for j in 1:length(F0_vec)]
    vdP_stream_forced(u) = Point2f( constant_forced_vdP(u, p, 0.0) )

    for i in 1:length(F0_vec)
        F0 = F0_vec[i]
        p["F₀"] = F0
        sp = streamplot!(
            axs[i],
            vdP_stream_forced,
            -3 .. 3,
            -3 .. 3,
            gridsize = (32, 32),
        )
        xlims!(axs[i], (-3, 3))
        ylims!(axs[i], (-3, 3))
    end
elseif case == 5
    fig = Figure(resolution = (1500, 500), font = srcdir("cmunrm.ttf"), fontsize = 20)
    F0_vec = 1:.1:1.2
    axs = [Axis( fig[1,j], xlabel = L"$x_{1}$ [m]", ylabel = L"$x_{2}$ [m/s]", ) for j in 1:length(F0_vec)]
    vdP_stream_forced(u) = Point2f( constant_forced_vdP(u, p, 0.0) )

    for i in 1:length(F0_vec)
        F0 = F0_vec[i]
        p["F₀"] = F0
        sp = streamplot!(
            axs[i],
            vdP_stream_forced,
            0.5 .. 1.5,
            -0.5 .. 0.5,
            gridsize = (32, 32),
        )
        xlims!(axs[i], (0.5, 1.5))
        ylims!(axs[i], (-0.5, 0.5))
    end
end

xlims!(ax, (-3, 3))
ylims!(ax, (-3, 3))
# axislegend(ax)
save_fig( plotsdir("vdp/"), "stream_traj", "both", fig )