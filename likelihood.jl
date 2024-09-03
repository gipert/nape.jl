using Measurements
using TypedTables
using Distributions
using DensityInterface


Qββ = 2039.04 # keV
N_A = 6.022E23
m_76 = 75.6E-3 # kg/mol
ΔE = 240 # keV

# Tuple{Real, Real, Vector{Real}, Vector{Real}, Real}
# k-parameters are vectors indexed by event indices (timestamp ordering)
ModelParameters = NamedTuple{(:Γ12, :B, :Δk, :σk, :α)}

# TODO: check that all events are found in the partitions table

# experiment = data set with shared background index parameter B
function make_exp_likelihood(events::Table, partitions::Table)
    DensityInterface.logfuncdensity(
        p::ModelParameters -> begin
            logL = 0

            # loop over partitions
            for part in partitions
                # calculate some useful stuff
                # same notation as in LNote 24-006
                # FIXME: remove math with big/small numbers
                μsk = log(2) * N_A * part.exposure * (part.ϵk.val + p.α * part.ϵk.err) * p.Γ12 / m_76
                μbk = p.B * ΔE * part.exposure
                μk = μbk + μsk

                # TODO: Julia syntax to explode a Measurement in two arguments?
                logL += logpdf(Poisson(μk), length(E)) + logpdf(Normal(), p.α)

                # FIXME: event_idx ordering: order of partitions table
                for i in part.event_idxs
                    logL += - log(μk) +
                            + log(μbk/ΔE + μsk * pdf(Normal(Qββ + p.Δk[i], p.σk[i]), events[i].energy))
                end
            end

            return logL
        end
    )
end
