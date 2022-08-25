using .Threads
include("utils.jl")

mutable struct SlicedScatterStructs
    avec::Vector{Float64}
    Fvec::Vector{Float64}
end

function get_scatter!(node, plot_type, sss, solve_ode, grid_dict, sol_dict)
    println( string("Getting results for ", plot_type, " = $node") )

    Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

    if plot_type == "Δx"
        Δx = node
    elseif plot_type == "D"
        d = D2d(node, p)
        p["d"] = d
        Δx = p["Δx"]
    elseif plot_type == "μ"
        p["μ"] = node
        Δx = p["Δx"]
    elseif plot_type == "σ"
        p["σ"] = node
        Δx = p["Δx"]
    end
    x₀ = get_x₀(p, Δx)
    t_buffer = p["t_buffer"]
    method = "last"     # can also use "mean"

    @showprogress for a in sss.avec
        p["aF"] = a
        if p["fixed_tend"] == true
            local tend = p["F_crit"] / minimum(sss.avec) + t_buffer
        else
            local tend = p["F_crit"] / a + t_buffer
        end
        local tspan = (0.0, tend)   # Simulation time span.
        
        for Fmax in sss.Fvec

            p["Fmax"] = Fmax
            p["t₂"] = p["t₁"] + Fmax / a
            append!(Fmax_scat, Fmax)
            append!(a_scat, a)

            local sol = solve_ode(x₀, tspan, p)
            # sol_dict[string(Fmax, "  ", a)] = sol
            if method == "last"
                append!(x_scat, norm(last(sol), 2))
            else
                nend = 100
                end_sol = sol(range(tend - t_buffer / 10, stop = tend, length = nend))
                end_position = reverse([end_sol.u[end-i][1] for i = 0:nend-1])
                append!(x_scat, abs(mean(end_position)))
            end
        end
    end
    grid_dict[string(node)] = [Fmax_scat, a_scat, x_scat]
end

function append_fparams!(F, a, p)
    p["t₂"] = p["t₁"] + p["Fmax"] / p["aF"]
    append!(F, p["Fmax"])
    append!(a, p["aF"])
end

function threaded_stochastic_scatter(node, sss, solve_ode, nmembers)
    println("Getting results for σ = $node")
    p["σ"] = node
    Δx = 0.6
    x₀ = get_x₀(p, Δx)
    Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)

    # Fix end time to time required for lowest slope to reach Fcrit, plus a buffer time.
    tend = p["F_crit"] / minimum(sss.avec) + p["t_buffer"]
    tspan = (0.0, tend)

    @showprogress for a in sss.avec
        p["aF"] = a
        for Fmax in sss.Fvec
            p["Fmax"] = Fmax
            append_fparams!(Fmax_scat, a_scat, p)
            count_on_thread = init_threads(0.0)
            @threads for i in 1:nmembers
                local sol = solve_ode(x₀, tspan, p)
                add_to_count = ( norm(last(sol), 2) > p["xₜ"] )
                @inbounds count_on_thread[threadid()] += add_to_count
            end
            count = sum(count_on_thread)
            empirical_prob = count / nmembers
            append!(x_scat, empirical_prob)
        end
    end
    return [Fmax_scat, a_scat, x_scat]
end

function init_threads(A)
    return zeros(eltype(A), nthreads())
end

function update_R!(R, A)
    @threads for i in eachindex(A)
        @inbounds R[threadid()] += A[i]
    end
end

# @btime clean_stochastic_scatter(0.5, sss, solve_plo, 100)
# @btime threaded_stochastic_scatter(0.5, sss, solve_plo, 100)