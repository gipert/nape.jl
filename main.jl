using BAT
using Distributions
using IntervalSets: (..)
using Plots

include("gerda.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

events = read_events_gerdaII()
partitions = read_partitions_gerdaII()
likelihood = make_exp_likelihood(events, partitions)

events = add_partition_info(events, partitions)
partitions = add_event_idxs(partitions, events)

prior = distprod(
    Γ12 = 1E-28..1E-26,
    B = 0..1,
    Δk = [Normal(v.val, v.err) for v in events.Δk],
    σk = [truncated(Normal(v.val, v.err), lower=0) for v in events.σk],
    α = 0..10,
)

posterior = PosteriorMeasure(likelihood, prior)
samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

bat_report(samples)
