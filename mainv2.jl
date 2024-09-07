using BAT
using Distributions
using IntervalSets: (..)
using Plots

# using Debugger
# break_on(:error)

include("gerda.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

experiments = (:gerdaII, :legend200)

full_likelihood = logfuncdensity(
    p -> prod([likelihood(get_data(exp)..., (p.Γ12, p[exp]...)) for exp in experiments])
)

prior = distprod(
    # FIXME: why does lower=0 crash BAT?
    Γ12 = 0.01..1,
    gerdaII = distprod(
        B = 1E-5..1E-2, # cts / keV kg yr
        Δk = [Normal(v.val, v.err) for v in get_data(:gerdaII).events.Δk],
        # FIXME: BAT is bugged, crashes with truncated distributions
        σk = [Normal(v.val, v.err) for v in get_data(:gerdaII).events.σk],
        α = 0..1,
    ),
    legend200 = distprod(
        B = 1E-5..1E-2, # cts / keV kg yr
        Δk = [Normal(v.val, v.err) for v in get_data(:legend200).events.Δk],
        # FIXME: BAT is bugged, crashes with truncated distributions
        σk = [Normal(v.val, v.err) for v in get_data(:legend200).events.σk],
        α = 0..1,
    ),
)

posterior = PosteriorMeasure(full_likelihood, prior)

@time samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)
).result

bat_report(samples)
