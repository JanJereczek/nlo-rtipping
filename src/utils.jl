include("forcing.jl")

# Create a meshgrid with input vectors x and y.
function meshgrid(x, y)
    return (repeat(x, outer = length(y)), repeat(y, inner = length(x)))
end

# Compute sigmoid with s=steepness parameter, x_inf=inflection point, x=variable
function sigmoid(s::Float64, x_inf::Float64, x::Float64)
    1 / (1 + exp(-s * (x - x_inf)))
end

function normalise(x)
    # println(mean(x), "    ",var(x))
    return (x .- mean(x)) ./ sqrt(var(x))
end

function normalise_ext(x, μ, σ2)
    return (x .- μ) ./ sqrt(σ2)
end

# Convert angular frequency to frequency.
function ω2f(ω::Any)
    return ω ./ (2 * π)
end

# Convert frequency to angular frequency.
function f2ω(f::Any)
    return (2 * π) .* f
end

function get_x₀(p, Δx)
    return [p["xeq1"] - Δx, 0.0]
end

##############################################################################
############################ Fourier Transforms ##############################
##############################################################################
# Notice: some terms of the computed spectra are 0 for any f != 0. 
# (essentially dirac delta and its derivative).
# As log(0) = -∞, it won't be of importance and we truncate these terms.

# Compute Fourier transform of saturated ramp with parameters t1, t2, Fmax, a.
# Output image of input ω.
function ft_satramp(t1::Float64, t2::Float64, a::Any, ω::Any)
    s = im .* ω
    return a ./ s.^2 .* (exp.(-t1 .* s) .- exp.(-t2 .* s))
end

# Compute Fourier transform of step with amplitude Fstep.
function ft_step(Fstep::Float64, ω::Any)
    return Fstep ./ (im .* ω)
end

# Compute Fourier transform of saturated ramp + step.
function ft_stepramp(
    t1::Float64,
    t2::Float64,
    Fstep::Float64,
    a::Any,
    ω::Any,
)
    return ft_satramp(t1, t2, a, ω) .+ ft_step(Fstep, ω)
end

function get_spectral_power(get_resp1, fspace1, get_resp2)
    return get_resp1(fspace1) .* get_resp2(fspace1)
end

function save_fig(prefix, filename, extension, fig)
    if (extension == "pdf") || (extension == "both")
        save(string(prefix, filename, ".pdf"), fig)
    end
    if (extension == "png") || (extension == "both")
        save(string(prefix, filename, ".png"), fig)
    end
end



# function ft_satramp(t1::Float64, t2::Float64, Fmax::Float64, a::Any, ω::Any)
#     k = ω ./ a
#     return (Fmax  ./ abs.(a)) ./ (k .^ 2) .*
#            (exp.(-im .* k .* t2) - exp.(-im .* k .* t1))
# end