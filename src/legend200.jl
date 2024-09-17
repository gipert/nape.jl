using HDF5
using LegendHDF5IO
using LegendDataManagement
using TypedTables
using JSON
using Measurements
using IntervalSets: Interval, (..)

function read_events_legend200()::Table
    cfg = LegendDataConfig("data/legend200/config.json")
    ldata = LegendData(cfg.setups.l200, :l200)
    lmeta = ldata.metadata

    data = lh5open("data/legend200/l200-nu24-phy-tier_skm.lh5")["skm"][:]

    E = data.geds.energy
    # cannot code a union (∪) of disjoint intervals with IntervalSets?
    e_sel = E .∈ [1930..2099] .|| E .∈ [2109..2114] .|| E .∈ [2124..2190]

    data = data[
        e_sel .&&
        .!data.coincident.muon_offline .&&
        .!data.coincident.spms .&&
        data.geds.psd.is_bb_like
    ]

    detector = [
        # thanks florian
        channelinfo(ldata, FileKey("l200-p03-r000-phy-20230312T043356Z"), Int(rawid)).detector.label
        for rawid in data.geds.rawid
    ]

    # data types should be correct
    return Table(
        timestamp=round.(Int64, data.trigger.timestamp),
        detector=detector,
        energy=data.geds.energy
    )
end

function read_partitions_legend200()::Table
    span, detector, exposure, ϵk, Δk, σk = _make_partitions_columns()

    cfg = LegendDataConfig("data/legend200/config.json")
    lmeta = LegendData(cfg.setups.l200, :l200).metadata

    detectors = lmeta.datasets.ovbb_partitions_pars
    # loop over detectors
    for (hpge, detdata) in detectors
        # apply defaults
        hpge == :default && continue
        _default = :default in keys(detectors) ? detectors.default : Dict()
        detdata = merge(_default, detdata)

        for (pname, pardata) in detdata
            # apply defaults
            pname == :default && continue
            _default = :default in keys(detdata) ? detdata.default : Dict()
            pardata = merge(_default, pardata)

            push!(span, Interval(pardata.span_in_utc_s...))

            push!(detector, hpge)

            exp = pardata.livetime_in_s / 31536000 *
                  lmeta.hardware.detectors.germanium.diodes[hpge].production.mass_in_g / 1000
            push!(exposure, exp)

            enrichment = lmeta.hardware.detectors.germanium.diodes[hpge].production.enrichment
            eff = prod([v.val for v in values(pardata.ovbb_acceptance)]) * enrichment.val
            # FIXME: I'm sure there is a Julia builtin for sumsquared
            σ_eff = √(sum([v.unc for v in values(pardata.ovbb_acceptance)].^2) + enrichment.unc^2)
            push!(ϵk, eff ± σ_eff)

            # NOTE: seems like the definition of bias in the config is opposite to ours (the likelihood)
            push!(Δk, - pardata.energy_bias_in_keV.val ± pardata.energy_bias_in_keV.unc)
            push!(σk, (pardata.fwhm_in_keV.val ± pardata.fwhm_in_keV.unc) / 2.355)
        end
    end

    return Table(span=span, detector=detector, exposure=exposure, ϵk=ϵk, Δk=Δk, σk=σk)
end
