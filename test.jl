using DensityInterface
using Debugger

include("gerda.jl")
include("legend200.jl")
include("likelihood.jl")
include("tools.jl")

events, partitions = get_data(:gerdaII)

@time @show likelihood(
    events,
    partitions,
    (Γ12=0.5, B=0.0004, Δk=fill(0.01, length(events)), σk=fill(2.5, length(events)), α=0.01)
)

@show likelihood(
    events,
    partitions,
    (Γ12=0, B=0.0004, Δk=fill(0.01, length(events)), σk=fill(2.5, length(events)), α=0.01)
)
