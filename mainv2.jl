using BAT
using Distributions
using IntervalSets: (..)
using Plots

# using Debugger
# break_on(:error)

include("gerda.jl")
include("legend200.jl")
include("majorana.jl")
include("likelihood.jl")
include("tools.jl")

experiments = (:gerdaII, :legend200)

function getpars(pars::NamedTuple, experiment::Symbol)::NamedTuple
    sel = filter(x -> startswith(string(x.first), "$(experiment)_"), pairs(pars))
    return NamedTuple(Symbol(chopprefix(string(k), "$(experiment)_")) => v for (k, v) in pairs(sel))
end

full_likelihood = logfuncdensity(
    # non-global parameter names are prefixed by the experiment id, need to
    # strip it off to make the likelihood work
    p -> prod([likelihood(get_data(exp)..., (Γ12=p.Γ12, getpars(p, exp)...)) for exp in experiments])
)

# FIXME: BAT is bugged, crashes with truncated distributions
prior = distprod(
    Γ12 = 0.01..1,
    gerdaII_B = 1E-5..1E-2, # cts / keV kg yr
    gerdaII_Δk = [Normal(v.val, v.err) for v in get_data(:gerdaII).events.Δk],
    gerdaII_σk = [Normal(v.val, v.err) for v in get_data(:gerdaII).events.σk],
    gerdaII_α = 0..1,
    legend200_B = 1E-5..1E-2, # cts / keV kg yr
    legend200_Δk = [Normal(v.val, v.err) for v in get_data(:legend200).events.Δk],
    legend200_σk = [Normal(v.val, v.err) for v in get_data(:legend200).events.σk],
    legend200_α = 0..1,
)

posterior = PosteriorMeasure(full_likelihood, prior)

@time samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)
).result

bat_report(samples)
