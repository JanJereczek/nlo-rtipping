

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