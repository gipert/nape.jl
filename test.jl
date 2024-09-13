using DensityInterface
using BenchmarkTools
using Debugger

include("gerda.jl")
include("majorana.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

experiments = (
    :gerdaI_golden, :gerdaI_silver, :gerdaI_bege, :gerdaI_extra,
    :gerdaII,
    :majorana_DS0, :majorana_mod1, :majorana_mod2,
    :legend200
)

pars = (
    Γ12=0.5,

    gerdaII_B=4E-4,
    gerdaII_Δk=fill(0.01, 13),
    gerdaII_σk=fill(2.5, 13),
    gerdaII_α=0.01,

    legend200_B=2E-4,
    legend200_Δk=fill(0.01, 7),
    legend200_σk=fill(2.4, 7),
    legend200_α=0.02
)

data = Dict(exp => get_data(exp) for exp in experiments)

for _logl in (loglikelihood, loglikelihood_experimental)

    @info "Benchmarking: " _logl

    display(@benchmark $_logl(data[:gerdaII]..., pars))
    @show _logl(data[:gerdaII]..., pars)
end

# @show getpars(pars, :gerdaII)
# @show getpars(pars, :legend200)

# full_loglikelihood = logfuncdensity(
#     p -> sum([loglikelihood(data[exp]..., (Γ12=p.Γ12, getpars(p, exp)...)) for exp in experiments])
# )

# display(@benchmark logdensityof(full_loglikelihood)(pars))
# @show logdensityof(full_loglikelihood)(pars)
