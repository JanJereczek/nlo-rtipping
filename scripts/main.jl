# Import packages.
using DifferentialEquations, CairoMakie, Colors, NFFT, Interpolations, ProgressMeter;

# Include self-written scripts.
include("forcing.jl")
include("mechanics.jl")
include("nonlinear_oscillator.jl")
include("parameters.jl")
include("plots.jl")
include("utils.jl")

# Define saving options.
save_prefixes = ["plots_v0.2/tmp/", "plots_v0.2/pertx0/", "plots_v0.2/equilx0/"]
prefix = save_prefixes[1]

# Define plotting options.
plot_characteristics_bool = true
plot_bifurcation_bool = true
plot_stream_bool = true
plot_tip_grid_bool = false
plot_response_bool = true

# Define further options.
compute_ft_analytically = true
noisy_forcing = false

# Parameter dictionnary.
p = load_parameters()
p["F_type"] = "constant"
Δx = 0.44
p["t₁"] = 0.0
p["xeq1"] = 1.0

# Control Run.
if prefix == "plots/equilx0/"
    x₀ = [p["xeq1"], 0.0]               # Start at equilibrium.
    plot_equil_control(p, prefix) # Plot control run of x0 = 0.99*xtip.
else
    x₀ = [p["xeq1"] - Δx, 0.0]               # Start at Δx=1m from equilibrium.
end

# Get characteristics of spring n°2.
x₁_vec = range(0.0, stop = 2.0, length = 1000)
c₂ = get_nl_stiffness.(x₁_vec)

# Get frequency characterstics.
f_res, ω_res = round.(get_resonance_freq(p); digits = 3)

ω_vec = 10 .^ (range(-3.0, stop = 2.0, length = 1001))
f_vec = ω_vec ./ (2 * π)
amp_resp = get_bode_amp.(ω_vec)
get_ampresp = LinearInterpolation(f_vec, amp_resp)

ix_res = argmax(amp_resp)
fres_vec = f_vec[amp_resp.>1.1]
ampres_vec = amp_resp[amp_resp.>1.1]
f1_res, f2_res = round(fres_vec[1]; digits = 3), round(fres_vec[end]; digits = 3)
println("Resonance freq (Hz) - lower bound, peak and upper bound: ", f1_res, "  ", f_res, "  ", f2_res)

# Plot the characteristics
if plot_characteristics_bool
    plot_characteristics(x₁_vec, c₂, f_vec, amp_resp, prefix)
end

# Here we might compute the spectral power in future versions.
# spectral_power = get_spectral_power(get_ampresp, f_vec, get_ampresp, f_vec)
# println(spectral_power)

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

# Choose forcing type.
if noisy_forcing
    p["F_type"] = "noisy"
    p["σ"] = 1.0
    prefix = string(prefix, "noisy_")
else
    p["F_type"] = "ramp"
end

# Sampled values of the tipping grid.
nF, na = 50, 50
F_llim, F_ulim = -15, 2     # linscale with ref = F_crit
a_llim, a_ulim = -2, 3      # logscale with ref = 10^0
Fvec = round.( p["F_crit"] .+ range(F_llim, stop = F_ulim, length = nF); digits=5)
avec = round.( 10 .^ (range(a_llim, stop = a_ulim, length = na)); digits=5)
n_int = 100
cb_maps = [ :rainbow, :viridis ]
cb_limits = [(1.5, 3.1), (50.0, 100.0)]

# Initialise figure of tipping grid.
Fmax_scat, a_scat, resfreq_scat, x_scat = zeros(0), zeros(0), zeros(0), zeros(0)
grid_fig = Figure(resolution = (1600, 800), fontsize = 18)
resp_fig = Figure(resolution = (1500, 1000), fontsize = 18)

##################################################################
############### Tipping grid and system response #################
##################################################################

if plot_tip_grid_bool | plot_response_bool

    grid_axs = init_grid_axs(grid_fig)
    resp_axs = init_response_axs(resp_fig, f1_res, f2_res)

    for a in avec
        p["aF"] = a
        local tend = p["F_crit"] / a + 10
        local tspan = (0.0, tend)   # Simulation time span.
        local t_vec = range(tspan[1], stop = tspan[2], step = 0.1)
        # msg = string("Getting the tipping grid for a=", a, "...")
        # @showprogress msg 
        println(a)
        for Fmax in Fvec

            p["F_max"] = Fmax
            p["t₂"] = p["t₁"] + Fmax / a
            append!(Fmax_scat, Fmax)
            append!(a_scat, a)

            if compute_ft_analytically
                Ft = ft_pert_ramp(p["t₁"], p["t₂"], Fmax, a, ω_vec, - get_c(p["xeq1"]-Δx)*Δx )
            else
                Ft = fft_ramp(t_vec, plan)
            end
            
            if plot_tip_grid_bool
                local sol = solve_nlo(x₀, tspan, p)
                append!(x_scat, last(sol)[1])

                get_Ft = LinearInterpolation(f_vec, abs.(Ft))
                # A_f_res_lininterp = get_Ft(f_res)
                # A_f_res_argmax = abs.(Ft)[ix_res]

                fres_sample = range(f1_res, stop=f2_res, length=n_int)
                freq_product = get_Ft.( fres_sample ) .* get_ampresp.( fres_sample )
                A_f_res_integrated = sum( freq_product ) / n_int
                append!(resfreq_scat, A_f_res_integrated)
            end

            if plot_response_bool
                f_range = 45 .. 45.2
                a_range = 7 .. 7e1
                if (Fmax ∈ f_range) & (a ∈ a_range)
                    show_response(x₀, Fmax, a, t_vec, p, f_vec, abs.(Ft), resp_axs)
                end
            end

        end
    end

    if plot_tip_grid_bool

        hm1 = scatter!(grid_axs[1], a_scat, Fmax_scat, color = x_scat, colormap = cb_maps[1], colorrange = cb_limits[1])
        Colorbar(grid_fig[1, 2], hm1)

        hm2 = scatter!(grid_axs[2], a_scat, Fmax_scat, color = resfreq_scat, colormap = cb_maps[2], colorrange = cb_limits[2])
        Colorbar(grid_fig[1, 4], hm2)

        save_fig(prefix, "tipping_grid", "both", grid_fig)
    end

    if plot_response_bool
        resp_fig[1:2, 3] = Legend(resp_fig, resp_axs[4], framevisible = false)
        save_fig(prefix, "system_response", "both", resp_fig)
    end

end
