using BAT
using Distributions
using IntervalSets: (..)
using Optim
using Printf

include("data.jl")
include("legend200.jl")
include("majorana.jl")
include("likelihood.jl")
include("tools.jl")

# :gerdaI_golden, :gerdaI_silver, :gerdaI_bege, :gerdaI_extra,
# :gerdaII,
# :majorana_DS0, :majorana_mod1, :majorana_mod2,
# :legend200
data = get_data(:gerdaII)

# let block because accessing global stuff is slow?
full_loglikelihood = let data=data
    logfuncdensity(p -> loglikelihood_1bkg(data..., p...))
end

prior = distprod(
    Γ12 = 0..2,
    B = 1E-5..1E-2, # cts / keV kg yr
    Δk = [Normal(v.val, v.err) for v in data.events.Δk],
    σk = [truncated(Normal(v.val, v.err), lower=0) for v in data.events.σk],
    α = truncated(Normal(), lower=-2),
)

posterior = PosteriorMeasure(full_loglikelihood, prior)

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
@printf "T_12 > %.3g (90%% CI)" 1E26 ./ quantile(samples, 0.9).Γ12
