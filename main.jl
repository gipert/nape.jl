using BAT
using Optim
using Distributions
using IntervalSets: (..)
using Printf
using Plots

# using Debugger
# break_on(:error)

include("gerda.jl")
include("legend200.jl")
include("majorana.jl")
include("likelihood.jl")
include("tools.jl")

experiments = (
    :gerdaI_golden, :gerdaI_silver, :gerdaI_bege, :gerdaI_extra,
    :gerdaII,
    :majorana_DS0, :majorana_mod1, :majorana_mod2,
    :legend200
)

# cache data in memory
data = Dict(exp => get_data(exp) for exp in experiments)

# NOTE: it would be nice to further nest the parameters structure and segregate
# the non global parameters to a dedicated named tuple (one per experiment).
# this is unfortunately not possible yet with BAT, so I am going to just prefix
# the parameter names with the experiment name

# build the combined loglikelihood by summing the partial loglikelihoods
loglikelihood = let data=data, experiments=experiments, logl=loglikelihood_1bkg, _get=getpars
    DensityInterface.logfuncdensity(
        # non-global parameter names are prefixed by the experiment id, need to
        # strip it off to make the likelihood work -> see getpars()
        p -> sum(
            [
                logl(data[exp]..., _get(p, exp, (:Γ12, :B, :Δk, :σk, :α))...)
                for exp in experiments
            ]
        )
    )
end

ϵk = (exp, k) -> getfield.(data[exp].partitions.ϵk, k)

# get minimum value that α can assume
α_min = minimum([minimum(ϵk(exp, :val) ./ ϵk(exp, :err)) for exp in experiments])

# automatize building prior distributions for an experiment
make_exp_priors = experiment -> begin
    # get the energy scale/resolution partition data to construct the priors
    # to do this, we get the partitions for which there is at least one observed event
    # the ordering is defined by the events table
    _d = data[experiment]
    _pdata = _d.partitions[unique(_d.events.part_idx)]

    # prefix experiment label to parameter
    keys = (Symbol("$(experiment)_$k") for k in (:B, :Δk, :σk))

    values = (
#= B =# 1E-5..1E-1, # cts / keV kg yr
#= Δ =# [Normal(v.val, v.err) for v in _pdata.Δk],
        # energy resolution cannot be negative!
#= σ =# [truncated(Normal(v.val, v.err), lower=0) for v in _pdata.σk],
    )

    # zip
    (; (keys .=> values)...)
end

# prior & posterior definition

prior = distprod(;
    # the 0vbb half-life is a global parameter
    Γ12 = 0..5, # 10^-26 yr^-1
    α = truncated(Normal(), lower=-α_min),
    # all the other parameters are experiment specific
    merge([make_exp_priors(exp) for exp in experiments]...)...
)

posterior = PosteriorMeasure(loglikelihood, prior)

@info "starting to generate posterior samples"

@time samples = bat_sample(
    posterior,
    MCMCSampling(
        mcalg=MetropolisHastings(),
        nsteps=100_000, nchains=4,
        burnin=MCMCMultiCycleBurnin(
            nsteps_per_cycle=10000,
            max_ncycles=30
        )
    )
).result

# refined global mode search
globalmode = bat_findmode(
    posterior,
    OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([mode(samples)]))
).result

# console printing

println(bat_report(samples))

for exp in experiments
    sym = Symbol("$(exp)_B")
    B_68 = BAT.smallest_credible_intervals(getproperty(samples.v, sym))[1] / 1E-4
    @printf "[%s] BI = %.3g [%.3g, %.3g] (68%% CI) 10^-4 cts / keV kg yr\n" exp globalmode[sym] / 1E-4 B_68.left B_68.right
end

@printf "T_12 > %.3g (90%% CI)" 1E26 ./ quantile(samples, 0.9).Γ12
