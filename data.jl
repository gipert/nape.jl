using HDF5
using LegendHDF5IO
using TypedTables
using JSON
using Measurements

# events: Table with columns (timestamp, detector, energy). timestamp and
# detector must match the format in the partitions table
#
# partitions: Table with columns (
#     span::Tuple{Real},
#     detector::String,
#     exposure::Real, in kg yr
#     ϵk::Measurement,
#     Δk::Measurement, in keV
#     σk::Measurement, in keV
# )

function read_events_l200()::Table
    data = lh5open("data/legend200/l200-nu24-phy-tier_skm.lh5")["skm"][:]

    E = data.geds.energy
    # todo: Julia syntax to check if in interval? IntervalSets does not seem to support .∈
    e_sel = (E .> 1930 .&& E .< 2099) .|| (E .> 2109 .&& E .< 2114) .|| (E .> 2124 .&& E .< 2190)

    data = data[e_sel .&& .!data.coincident.muon_offline .&& .!data.coincident.spms .&& data.geds.psd.is_bb_like]

    return Table(
        timestamp=data.trigger.timestamp,
        # detector=data.geds.rawid,
        # FIXME: get from metadata
        detector=["V09374A", "V04199A", "B00061A", "P00661C", "V00048A", "V08682B", "V07647A"],
        energy=data.geds.energy
    )
end

# IN PROGRESS
function read_partitions_l200()::Table
    span = []; detector = []; exposure = []
    ϵk = []; Δk = []; σk = []

    # FIXME: switch to yaml
    detectors = JSON.parsefile("data/legend200/ovbb_partitions_pars.json")
    for (hpge, detdata) in detectors
        hpge == "default" && continue
        detdata = merge(get(detectors, "default", Dict()), detdata)

        for (pname, pardata) in detdata
            pname == "default" && continue
            pardata = merge(get(detdata, "default", Dict()), pardata)

            push!(span, Tuple(pardata["span_in_utc_s"]))
            push!(detector, hpge)
            # exp = pardata["livetime_in_s"] / 31536000 *
            #       * lmeta.hardware.detectors.germanium.diodes[hpge].production.mass_in_g / 1000
            push!(exposure, 0)

            # missing lmeta.hardware.detectors.germanium.diodes[hpge].production.enrichment
            eff = prod([v["val"] for v in values(pardata["ovbb_acceptance"])])
            # FIXME: I'm sure there is a Julia builtin for sumsquared
            σ_eff = √(sum([v["unc"] for v in values(pardata["ovbb_acceptance"])].^2))
            push!(ϵk, eff ± σ_eff)

            push!(Δk, pardata["energy_bias_in_keV"]["val"] ± pardata["energy_bias_in_keV"]["unc"])
            push!(σk, pardata["fwhm_in_keV"]["val"] ± pardata["fwhm_in_keV"]["unc"])
        end
    end

    return Table(span=span, detector=detector, exposure=exposure, ϵk=ϵk, Δk=Δk, σk=σk)
end

# tier4->Scan(
#   "timestamp:energy",
#   "!isTP && !isBL && !isMuVetoed && ((energy > 1930 && energy < 2099) || (energy > 2109 && energy < 2114) || (energy > 2124 && energy < 2190)) && !isPSDVetoed && !isLArVetoed && (multiplicity == 1) && psdIsEval && (datasetID != 3
#   )
# lookup detector name from "Instance" and gerda-metadata/detector-data
function read_events_gerdaII()::Table
    return Table(
        timestamp=[1455109448, 1457847659, 1472522222, 1475981084, 1480290460,
                   1485848926, 1503578885, 1509498133, 1516142805, 1533092526,
                   1539047354, 1566823934, 1568276649],
        detector=["ANG4", "GD61C", "GD35B", "ANG1", "GD35B", "GD91A", "GD76C",
                  "ANG1", "RG1", "GD61C", "IC74A", "ANG4", "GD32D"],
        energy=[1995.2452, 1958.6807, 2018.1346, 1950.9419, 2067.9735,
                2056.4280, 2042.0641, 1962.7372, 1957.5059, 1970.1398,
                2058.8776, 2015.8751, 2012.0643],
    )
end

function read_partitions_gerdaII()::Table
    span = []; detector = []; exposure = []
    ϵk = []; Δk = []; σk = []

    for block in values(JSON.parsefile("data/gerda/0vbb-analysis-parameters.json"))
        for (k, v) in block
            k in ["start_ts", "end_ts", "start_utc", "end_utc"] && continue

            push!(span, (block["start_ts"], block["end_ts"]))
            push!(detector, k)
            push!(exposure, v["exposure"])
            push!(ϵk, v["eff_tot"] ± v["eff_tot_sigma"])
            push!(Δk, v["bias"] ± v["energy_sigma"])
            push!(σk, (v["fwhm"] ± v["fwhm_sigma"])/2.35) # TODO: check number
        end
    end

    return Table(span=span, detector=detector, exposure=exposure, ϵk=ϵk, Δk=Δk, σk=σk)
end
