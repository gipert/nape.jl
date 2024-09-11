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

experiments = (:gerdaI_golden, :gerdaI_silver, :gerdaI_bege, :gerdaI_extra, :gerdaII)

# cache data in memory
data = Dict(exp => get_data(exp) for exp in experiments)

# NOTE: it would be nice to further nest the parameters structure and segregate
# the non global parameters to a dedicated named tuple (one per experiment).
# this is unfortunately not possible yet with BAT, so I am going to just prefix
# the parameter names with the experiment name

# TODO: getpars() is introducing latency...
full_loglikelihood = logfuncdensity(
    # non-global parameter names are prefixed by the experiment id, need to
    # strip it off to make the likelihood work -> see getpars()
    p -> sum([loglikelihood(data[exp]..., (Γ12=p.Γ12, getpars(p, exp)...)) for exp in experiments])
)

# automatize building prior distributions for an experiment
make_exp_priors = experiment -> begin
    # prefix experiment label to parameter
    keys = (Symbol("$(experiment)_$k") for k in (:B, :Δk, :σk, :α))
    values = (
#= B =# 1E-5..1E-2, # cts / keV kg yr
#= Δ =# [Normal(v.val, v.err) for v in data[experiment].events.Δk],
        # energy resolution cannot be negative!
#= σ =# [truncated(Normal(v.val, v.err), lower=0) for v in data[experiment].events.σk],
#= α =# -1..1,
    )

    # zip
    (; (keys .=> values)...)
end

# prior & posterior definition

prior = distprod(;
    # the 0vbb half-life is a global parameter
    Γ12 = 0..2, # 10^-26 yr^-1
    # all the other parameters are experiment specific
    merge([make_exp_priors(exp) for exp in experiments]...)...
)

posterior = PosteriorMeasure(full_loglikelihood, prior)

@info "starting to generate posterior samples"

@time samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=100_000, nchains=4)
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

@printf "T_12 > %.3g (90%% CI)" quantile(1E26 ./ samples.v.Γ12, 0.1)
