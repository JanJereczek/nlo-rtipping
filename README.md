# Rate-induced tipping of nonlinear oscillators

## Getting started

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> nlo-rtipping

It is authored by Jan Swierczek-Jereczek.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Published paper

This repository contains the code of [Swierczek-Jereczek et al. (2023), *Time-scale synchronisation of oscillatory responses can lead to non-monotonous R-tipping*, Scientific reports.](https://www.nature.com/articles/s41598-023-28771-1)

The paper consists of a case study showing that rate-induced tipping might not be monotonous with respect to the forcing rate, a somewhat counter-intuitive result. It furthermore shows that for oscillatory systems, this might happen due to the synchronisation of oscillatory response.