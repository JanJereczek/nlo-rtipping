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


# function clean_stochastic_scatter(node, sss, solve_ode, nmembers)
#     println("Getting results for σ = $node")

#     p["σ"] = node
#     Δx = 0.6
#     x₀ = get_x₀(p, Δx)

#     Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)
#     sol_dict = Dict()

#     # Fix end time to time required for lowest slope to reach Fcrit, plus a buffer time.
#     t_buffer = 100
#     tend = p["F_crit"] / minimum(sss.avec) + t_buffer
#     tspan = (0.0, tend)

#     @showprogress for a in sss.avec
#         p["aF"] = a
#         for Fmax in sss.Fvec
#             p["Fmax"] = Fmax
#             append_fparams!(Fmax_scat, a_scat, p)
#             count = 0.0f0
#             for i in 1:nmembers
#                 local sol = solve_ode(x₀, tspan, p)
#                 add_to_count = Float32( ( norm(last(sol), 2) > p["xₜ"] ) )
#                 count += add_to_count
#             end
#             empirical_prob = count / nmembers
#             append!(x_scat, empirical_prob)
#         end
#     end
#     return [Fmax_scat, a_scat, x_scat]
# end



# function stochastic_scatter(node, sss, solve_ode, nmembers)
#     Fmax_scat, a_scat, x_scat = zeros(0), zeros(0), zeros(0)
#     p["σ"] = node
#     Δx = 0.6
#     x₀ = get_x₀(p, Δx)
#     println("Getting results for σ = $node")
#     sol_dict = Dict()
#     t_buffer = 100
#     nend = 100
#     method = "last"

#     @showprogress for a in sss.avec
#         p["aF"] = a
#         local tend = p["F_crit"] / minimum(sss.avec) + t_buffer
#         local tspan = (0.0, tend)   # Simulation time span.

#         for Fmax in sss.Fvec

#             p["Fmax"] = Fmax
#             p["t₂"] = p["t₁"] + Fmax / a
#             append!(Fmax_scat, Fmax)
#             append!(a_scat, a)

#             count = 0.
#             for i in 1:nmembers
#                 local sol = solve_ode(x₀, tspan, p)
#                 if method == "last"
#                     add_to_count = Float32( ( norm(last(sol), 2) > p["xₜ"] ) )
#                 else
#                     end_sol = sol(range(tend - t_buffer / 10, stop = tend, length = nend))
#                     end_position = reverse([end_sol.u[end-i][1] for i = 0:nend-1])
#                     add_to_count = Float32( ( abs( mean(end_position) ) > p["xₜ"] ) )
#                 end
#                 count += add_to_count
#             end
#             empirical_prob = count / nmembers
#             append!(x_scat, empirical_prob)
#         end
#     end
#     return [Fmax_scat, a_scat, x_scat]
# end