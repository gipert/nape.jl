using Measurements
using TypedTables
using Distributions
using DensityInterface


Qββ = 2039.04 # keV
N_A = 6.022
m_76 = 75.6 # g/mol
ΔE = 240 # keV

# Tuple{Real, Real, Vector{Real}, Vector{Real}, Real}
# k-parameters are vectors indexed by event indices (timestamp ordering)
ModelParameters = NamedTuple{(:Γ12, :B, :Δk, :σk, :α)}

# experiment = data set with shared background index parameter B

# no for loops, only array programming
function likelihood(events::Table, partitions::Table, p::ModelParameters)::Float64
    P = partitions
    # same notation as in LNote 24-006
    # Γ12 will be in units of 10^-26 yr^-1
    μsk = p.Γ12 * log(2) * N_A * P.exposure .*
          (getfield.(P.ϵk, :val) .+ p.α * getfield.(P.ϵk, :err)) / m_76
    μbk = p.B * ΔE * P.exposure
    μk = μbk .+ μsk

    pid = events.part_idx

    return (
        # first part of the likelihood that does not depend on the observed events
        sum(logpdf.(Poisson.(μk), length.(P.event_idxs)) .+ logpdf(Normal(), p.α)) +
        # part that depends on the observed events
        # FIXME: remove abs() once problem with truncated distributions is fixed
        sum(
            -log.(μk[pid]) .+
            log.(μbk[pid]/ΔE .+ μsk[pid] .* pdf.(Normal.(Qββ .+ p.Δk, abs.(p.σk)), events.energy)))
    )

end

# slower version with for loops
function slow_likelihood(events::Table, partitions::Table, p::ModelParameters)::Float64
    logL = 0

    # loop over partitions
    for part in partitions
        # calculate some useful stuff
        # same notation as in LNote 24-006
        # Γ12 will be in units of 10^-26 yr^-1
        μsk = p.Γ12 * log(2) * N_A * part.exposure * (part.ϵk.val + p.α * part.ϵk.err) / m_76
        μbk = p.B * ΔE * part.exposure
        μk = μbk + μsk

        logL += logpdf(Poisson(μk), length(part.event_idxs)) + logpdf(Normal(), p.α)

        # event_idx ordering: order as in partitions table
        for i in part.event_idxs
            logL += - log(μk) +
                    + log(μbk/ΔE + μsk * pdf(Normal(Qββ + p.Δk[i], p.σk[i]), events[i].energy))
        end
    end

    return logL
end
