include(srcdir("utils.jl"))
include(srcdir("parameters.jl"))

p = load_parameters()

t₁ = 0.0
t₂ = 1e-10
Fmax = 100.0
a = Fmax / (t₂ - t₁)
ω = 10 .^ range(-3, stop=3, length=1000)

ft1 = ft_satramp(t₁, t₂, Fmax, a, ω)
ft2 = ft_step(Fmax, ω)

println( isapprox( sum( abs.(ft2 - ft1) ), 0; atol=1e-3) )

#