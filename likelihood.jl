using Measurements
using TypedTables
using Distributions
using DensityInterface

# events: Table with columns (timestamp, detector, energy). timestamp and
# detector must match the format in the partitions table
#
# partitions: Table with columns (
#     span::Interval,
#     detector::String,
#     exposure::Real, in kg yr
#     ϵk::Measurement,
#     Δk::Measurement, in keV
#     σk::Measurement, in keV
# )

Qββ = 2039.04 # keV
N_A = 6.02214
m_76 = 75.9214 # g/mol
ΔE = 240 # keV

# k-parameters are vectors indexed by event indices (timestamp ordering)
ModelParameters = NamedTuple{(:Γ12, :B, :Δk, :σk, :α)}
# Tuple{Real, Real, Vector{Real}, Vector{Real}, Real}

# experiment = data set with shared background index parameter B

# no for loops, only array programming
# NOTE: @inbounds, @views and @fastmath make a marginal difference
# NOTE: broadcasted operations should be already vectorized?
# NOTE: not hard-typing "p", otherwise problems when tuple fields order is different
function loglikelihood_experimental(events::Table, partitions::Table, p::NamedTuple)::Float64
    P = partitions
    # same notation as in LNote 24-006
    # Γ12 will be in units of 10^-26 yr^-1
    μsk = p.Γ12 * log(2) * N_A * P.exposure .*
          (getfield.(P.ϵk, :val) .+ p.α * getfield.(P.ϵk, :err)) / m_76
    μbk = p.B * ΔE * P.exposure
    μk = μbk .+ μsk

    # energy nuisance parameters are indexed according to index in events table
    pid = events.part_idx
    Δk = p.Δk[events.epar_idx]
    σk = p.σk[events.epar_idx]

    return (
        # first part of the likelihood that does not depend on the observed events
        sum(logpdf.(Poisson.(μk .+ eps.(μk)), length.(P.event_idxs))) +
        # part that depends on the observed events
        sum(
            -log.(μk[pid]) .+
            log.(μbk[pid]/ΔE .+ μsk[pid] .* pdf.(Normal.(Qββ .+ Δk, σk), events.energy))
        )
    )

end

# this function is needed to get the correct experiment parameters when the
# parameters in the prior structure are prefixed with the experiment name (when
# combining multiple experiments)
# TODO: speedup by building global lookup table
function getpars(pars::NamedTuple, experiment::Symbol)::NamedTuple
    cjosul = filter(x -> startswith(string(x.first), "$(experiment)_"), pairs(pars))
    return NamedTuple(Symbol(chopprefix(string(k), "$(experiment)_")) => v for (k, v) in pairs(cjosul))
end

# slower likelihood with for loops
function loglikelihood(events::Table, partitions::Table, p::NamedTuple)::Float64
    logL = 0

    # loop over partitions
    for part in partitions
        # calculate some useful stuff
        # same notation as in LNote 24-006
        # Γ12 will be in units of 10^-26 yr^-1
        μsk = p.Γ12 * log(2) * N_A * part.exposure * (part.ϵk.val + p.α * part.ϵk.err) / m_76
        μbk = p.B * ΔE * part.exposure
        μk = μbk + μsk

        logL += logpdf(Poisson(μk + eps(μk)), length(part.event_idxs))

        # event_idx ordering: order as in partitions table
        for ev in events[part.event_idxs]
            logL += - log(μk) +
                    + log(μbk/ΔE + μsk * pdf(Normal(Qββ + p.Δk[ev.epar_idx], p.σk[ev.epar_idx]), ev.energy))
        end
    end

    return logL
end
