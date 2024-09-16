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

data = Dict(exp => get_data(exp) for exp in experiments)

pars = (
    Γ12=0.5,
    α=0.01,

    gerdaI_golden_B=1E-2,
    gerdaI_golden_Δk=fill(-0.01, 46),
    gerdaI_golden_σk=fill(1.8, 46),

    gerdaI_silver_B=1E-3,
    gerdaI_silver_Δk=fill(0.01, 10),
    gerdaI_silver_σk=fill(1.8, 10),

    gerdaI_bege_B=5E-4,
    gerdaI_bege_Δk=fill(0.02, 3),
    gerdaI_bege_σk=fill(1.1, 3),

    gerdaI_extra_B=3E-4,
    gerdaI_extra_Δk=fill(-0.02, 2),
    gerdaI_extra_σk=fill(1.7, 2),

    gerdaII_B=4E-4,
    gerdaII_Δk=fill(0.01, 13),
    gerdaII_σk=fill(1.2, 13),

    majorana_DS0_B=2E-2,
    majorana_DS0_Δk=fill(0.03, 7),
    majorana_DS0_σk=fill(2.5, 7),

    majorana_mod1_B=1E-2,
    majorana_mod1_Δk=fill(-0.03, 75),
    majorana_mod1_σk=fill(1.1, 75),

    majorana_mod2_B=4E-2,
    majorana_mod2_Δk=fill(0.05, 14),
    majorana_mod2_σk=fill(1.1, 14),

    legend200_B=2E-4,
    legend200_Δk=fill(0.05, 7),
    legend200_σk=fill(1.1, 7),
)

loglikelihood = let data=data, _logl=loglikelihood_1bkg
    DensityInterface.logfuncdensity(
        p -> sum(
            [
                _logl(data[exp]..., getpars(p, exp, (:Γ12, :B, :Δk, :σk, :α))...)
                for exp in experiments
            ]
        )
    )
end

display(@benchmark logdensityof(loglikelihood)(pars))
@show logdensityof(loglikelihood)(pars)
