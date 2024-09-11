using DensityInterface
using BenchmarkTools
using Debugger

include("gerda.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

experiments = (:gerdaII, :legend200)

data = Dict(exp => get_data(exp) for exp in experiments)

display(
    @benchmark loglikelihood(data[:gerdaII]...,
        (Γ12=0.5, B=0.0004, Δk=fill(0.01, 13), σk=fill(2.5, 13), α=0.01)
    )
)
@show loglikelihood(data[:gerdaII]...,
    (Γ12=0.5, B=0.0004, Δk=fill(0.01, 13), σk=fill(2.5, 13), α=0.01)
)

pars = (
    Γ12=0.5,

    gerdaII_B=0.0004,
    gerdaII_Δk=fill(0.01, 13),
    gerdaII_σk=fill(2.5, 13),
    gerdaII_α=0.01,

    legend200_B=0.0004,
    legend200_Δk=fill(0.01, 7),
    legend200_σk=fill(2.5, 7),
    legend200_α=0.0
)

@show getpars(pars, :gerdaII)
@show getpars(pars, :legend200)

full_loglikelihood = logfuncdensity(
    p -> sum([loglikelihood(data[exp]..., (Γ12=p.Γ12, getpars(p, exp)...)) for exp in experiments])
)

display(@benchmark logdensityof(full_loglikelihood)(pars))
@show logdensityof(full_loglikelihood)(pars)
