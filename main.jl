using BAT
using Distributions
using IntervalSets: (..)
using Plots

include("data.jl")
include("likelihood.jl")

likelihood = make_exp_likelihood(read_events_gerdaII(), read_partitions_gerdaII())
n_events = length(read_events_gerdaII())

prior = distprod(
    Γ12 = 1E-28..1E-26,
    B = 0..1,
    Δk = fill(0..10, n_events),
    σk = fill(0..10, n_events),
    α = 0..10,
)

posterior = PosteriorMeasure(likelihood, prior)
samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4))
