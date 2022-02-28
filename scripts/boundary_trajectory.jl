using NLsolve, Latexify, DifferentialEquations, CairoMakie
include(srcdir("parameters.jl"))
include(srcdir("nonlinear_oscillator.jl"))

p = load_parameters()
w0 = p["ω₀1"]
D = p["D"]
wd = p["ωd"]
xₜ = p["xₜ"]
xmg = 1.0
A = xₜ - xmg

x₁(A, B, t) = exp(-D*w0*t) * ( A*cos(wd*t) + B*sin(wd*t) ) + xmg
x₂(A, B, t) = -D*w0*exp(-D*w0*t) * ( A*cos(wd*t) + B*sin(wd*t) ) + exp(-D*w0*t) * ( -A*wd*sin(wd*t) + B*wd*cos(wd*t) )
# fB(t) = A*( D*w0*cos(wd*t) + wd*sin(wd*t) ) / ( wd*cos(wd*t) - D*w0*sin(wd*t)  )
fB(t) = (exp(D*w0*t) * (xₜ - xmg) - A*cos(wd*t)) / (sin(wd*t))

function f!(F, x)
    F[1] = x₂(A, fB(x[1]), x[1])
end

initial_x = [2.0]
@time sol_nl = mcpsolve(f!, [1.0], [5.0], initial_x, reformulation = :smooth)

t_tip = sol_nl.zero[1]
# t_tip = 0.87
B = fB(t_tip)
println("B, t_tip = ", [B, t_tip])

x₀ = [x₁(A, B, 0), x₂(A, B, 0)]
println("Residual: ", [x₁(A, B, t_tip) - xₜ, x₂(A, B, t_tip)] )
p["F_type"] = "const"
p["F_const"] = 0
p["F_noise"] = true
p["σ"] = 0
p["dt"] = 0.01
sol = solve_nlo(x₀, (0, 5), p)
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, sol.t, sol[1, :])
vlines!(ax, [t_tip], color = :red)
fig