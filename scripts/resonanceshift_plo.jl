#################################################
#################################################
#################################################

using DrWatson
@quickactivate "nlo-rtipping"

# Import packages.
using CairoMakie, Colors, Interpolations, JLD2
using DifferentialEquations, LinearAlgebra, ProgressMeter

# Include self-written scripts.
include(srcdir("utils.jl"))
include(srcdir("integrate_tipgrid.jl"))
include(srcdir("plastic_oscillator.jl"))

##################

const α2 = 3e1
const c2_min = 10.
const c2_max = 190.
const g = 0.
const Δc2 = c2_max - c2_min
const x_inflection_spring2 = 1.

p = standard_parameters()
p["c₁"] = 10.0        # [kg/s²] stiffness of linear spring.
p["d"] = 5.0          # [kg/s] dampening factor. Criteria for non-oscillating behavior: d > 2*sqrt(m*c)
x₀_vec = [0.5, 1.5]
ω_vec = 10 .^ (-1:.01:1)

fig = Figure()
ax1 = Axis(fig[1,1], xscale = log10, yscale = log10)
for x₀ in x₀_vec
    local par(ω::Float64) = plastic_amp_response(ω, x₀)
    lines!(ax1, ω_vec, par.(ω_vec), label = L"$x_{0} = %$(x₀) \, \mathrm{m}$")
end
axislegend(ax1)

ω_crit_lofreq = get_ω₀( get_total_plastic_stiffness(1.5) )

t_support = 0:.01:30
Fharmonic = 2. * sin.( ω_crit_lofreq * t_support )
ax2 = Axis(fig[1,2])
ax4 = Axis(fig[2,2])

for x₀ in x₀_vec
    local F0 = po_equilibrium_F(x₀)
    local F = F0 .+ Fharmonic
    p["forcing"] = LinearInterpolation( t_support, F )
    local sol = solve_po([x₀, 0], extrema(t_support), p)
    local X = hcat(sol.u...)
    lines!(ax2, t_support, F)
    lines!(ax4, sol.t, X[1,:])
end
fig