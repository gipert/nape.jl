using BAT
using Distributions
using IntervalSets: (..)
using Printf

# using Debugger
# break_on(:error)

include("gerda.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

gerda = get_data(:legend200)

# let block because accessing global stuff is slow?
full_likelihood = let gerda=gerda
    logfuncdensity(p -> likelihood(gerda..., p))
end

prior = distprod(
    # FIXME: why does lower=0 crash BAT?
    Γ12 = 0.01..5,
    B = 1E-5..1E-2, # cts / keV kg yr
    Δk = [Normal(v.val, v.err) for v in gerda.events.Δk],
    # FIXME: BAT is bugged, crashes with truncated distributions
    σk = [Normal(v.val, v.err) for v in gerda.events.σk],
    α = 0..1,
)

posterior = PosteriorMeasure(full_likelihood, prior)

@time samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)
).result

println(bat_report(samples))

@printf "BI = %.3g (68%% CI)" mode(samples.v.B)
@printf "Result: T_12 > %.3g (90%% CI)" quantile(1E26 ./ samples.v.Γ12, 0.9)
