using BAT
using Distributions
using IntervalSets: (..)
using Printf

# using Debugger
# break_on(:error)

include("gerda.jl")
include("legend200.jl")
include("majorana.jl")
include("likelihood.jl")
include("tools.jl")

gerda = get_data(:gerdaII)

# let block because accessing global stuff is slow?
full_likelihood = let gerda=gerda
    logfuncdensity(p -> likelihood(gerda..., p))
end

prior = distprod(
    Γ12 = 0..5,
    B = 1E-5..1E-2, # cts / keV kg yr
    Δk = [Normal(v.val, v.err) for v in gerda.events.Δk],
    # FIXME: BAT is bugged, crashes with truncated distributions
    σk = [Normal(v.val, v.err) for v in gerda.events.σk],
    α = 0..1,
)

posterior = PosteriorMeasure(full_likelihood, prior)

@time samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=100_000, nchains=4)
).result

println(bat_report(samples))

globalmode = bat_findmode(
    posterior,
    OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([mode(samples)]))
).result

B_68 = BAT.smallest_credible_intervals(samples.v.B)[1] / 1E-4
@printf "BI = %.3g [%.3g, %.3g] (68%% CI) 10^-4 cts / keV kg yr\n" globalmode.B / 1E-4 B_68.left B_68.right
@printf "T_12 > %.3g (90%% CI)" quantile(1E26 ./ samples.v.Γ12, 0.1)
