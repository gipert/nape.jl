using TypedTables
using JSON
using Measurements
using IntervalSets: Interval, (..)

include("tools.jl")

# :gerdaI_golden
# :gerdaI_silver
# :gerdaI_bege
# :gerdaI_extra

function read_events_gerdaI(dataset::Symbol)::Table
    timestamp, detector, energy = _make_events_columns()

    lookup = Dict(
        :gerdaI_golden => "ph1_golden",
        :gerdaI_silver => "ph1_silver",
        :gerdaI_bege => "ph1_bege",
        :gerdaI_extra => "ph1_extra",
    )

    for line in readlines("data/gerda/I/events.txt")
        cols = split(line)
        if cols[end] == lookup[dataset]
            push!(timestamp, parse(Int, cols[2]))
            push!(detector, Symbol(cols[1]))
            push!(energy, parse(Float64, cols[3]))
        end
    end

    return Table(timestamp=timestamp, detector=detector, energy=energy)
end

function read_partitions_gerdaI(dataset::Symbol)::Table
    # https://doi.org/10.1103/PhysRevLett.111.122503
    if dataset == :gerdaI_golden
        return Table(span=[1320849782..1369143782], detector=[:coax],
                     exposure=[17.9], ϵk=[0.571 ± 0.033], Δk=[0.0 ± 0.2],
                     σk=[4.26 ± 0.14] ./ 2.355)
    elseif dataset == :gerdaI_silver
        return Table(span=[339945251..342589688], detector=[:coax],
                     exposure=[1.304], ϵk=[0.571 ± 0.033], Δk=[0.0 ± 0.2],
                     σk=[4.26 ± 0.14] ./ 2.355)
    elseif dataset == :gerdaI_bege
        return Table(span=[1320849782..1369143782], detector=[:bege],
                     exposure=[2.404], ϵk=[0.663 ± 0.022], Δk=[0.0 ± 0.2],
                     σk=[2.73 ± 0.19] ./ 2.355)
    elseif dataset == :gerdaI_extra
        return Table(span=[1370007782..1380548582], detector=[:coax],
                     exposure=[1.904], ϵk=[0.619 ± 0.07], Δk=[0.0 ± 0.2],
                     σk=[4.16 ± 0.19] ./ 2.355)
    end
end

# :gerdaII

# tier4->Scan(
#   "timestamp:energy",
#   "!isTP && !isBL && !isMuVetoed && ((energy > 1930 && energy < 2099) || (energy > 2109 && energy < 2114) || (energy > 2124 && energy < 2190)) && !isPSDVetoed && !isLArVetoed && (multiplicity == 1) && psdIsEval && (datasetID != 3)"
# )
# lookup detector name from "Instance" and gerda-metadata/detector-data
function read_events_gerdaII()::Table
    return Table(
        timestamp=Vector{Int64}([1455109448, 1457847659, 1472522222,
                                 1475981084, 1480290460, 1485848926,
                                 1503578885, 1509498133, 1516142805,
                                 1533092526, 1539047354, 1566823934,
                                 1568276649]),
        detector=Vector{Symbol}([:ANG4, :GD61C, :GD35B, :ANG1, :GD35B, :GD91A,
                                 :GD76C, :ANG1, :RG1, :GD61C, :IC74A, :ANG4,
                                 :GD32D]),
        energy=Vector{Float32}([1995.2452, 1958.6807, 2018.1346, 1950.9419,
                                2067.9735, 2056.4280, 2042.0641, 1962.7372,
                                1957.5059, 1970.1398, 2058.8776, 2015.8751,
                                2012.0643]),
    )
end

function read_partitions_gerdaII()::Table
    span, detector, exposure, ϵk, Δk, σk = _make_partitions_columns()

    for block in values(JSON.parsefile("data/gerda/II/0vbb-analysis-parameters.json"))
        for (k, v) in block
            k in ["start_ts", "end_ts", "start_utc", "end_utc"] && continue

            # partitions with null exposure? wtf?
            v["exposure"] <= 0 && continue

            push!(span, block["start_ts"]..block["end_ts"])
            push!(detector, Symbol(k))
            push!(exposure, v["exposure"])
            push!(ϵk, v["eff_tot"] ± v["eff_tot_sigma"])
            push!(Δk, v["bias"] ± v["energy_sigma"])
            push!(σk, (v["fwhm"] ± v["fwhm_sigma"])/2.355)
        end
    end

    return Table(span=span, detector=detector, exposure=exposure, ϵk=ϵk, Δk=Δk, σk=σk)
end
