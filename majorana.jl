using TypedTables
using Measurements
using IntervalSets: (..)

# :majorana_DS0
# :majorana_mod1
# :majorana_mod2

function read_events_majorana(dataset::Symbol)::Table
    timestamp = []; detector = []; energy = []

    lookup = Dict(:majorana_DS0 => "mjd-DS0", :majorana_mod1 => "mjd-mod1",:majorana_mod2 => "mjd-mod2")
    for line in readlines("data/majorana/events.txt")
        cols = split(line)
        if cols[end] == lookup[dataset]
            push!(timestamp, parse(Int, cols[2]))
            push!(detector, cols[1])
            push!(energy, parse(Float64, cols[3]))
        end
    end

    return Table(timestamp=timestamp, detector=detector, energy=energy)
end

function read_partitions_majorana()::Table
    # span = []; detector = []; exposure = []
    # ϵk = []; Δk = []; σk = []

    # int = s -> parse(Int, s)
    # flt = s -> parse(Float64, s)

    # open("data/majorana/parameter.txt") do f
    #     while true
    #         cols = split(readline(f))
    #         length(cols) == 0 && break

    #         push!(span, int(cols[2])..int(cols[3]))

    #         cols = split(readline(f))
    #         push!(detector, cols[2])
    #         push!(exposure, flt(cols[7]))
    #         push!(ϵk, flt(cols[5]) ± flt(cols[6]))
    #         push!(Δk, flt(cols[8]) ± flt(cols[9]))
    #         push!(σk, flt(cols[3]) ± flt(cols[4]))
    #     end
    # end

    # return Table(span=span, detector=detector, exposure=exposure, ϵk=ϵk, Δk=Δk, σk=σk)

    return Table(
        span=[1000..1050, 1100..1150, 1100..1150, 1100..1150],
        detector=["ppc1", "ppc1", "ppc2" ,"icpc"],
        exposure=[1.26, 49.6, 17.06, 3.122],
        ϵk=[00.611054 ± 0.023909, 0.62 ± 0.03, 0.62 ± 0.03, 0.592331 ± 0.033576],
        Δk=[0.0 ± 0.2, 0.0 ± 0.2, 0.0 ± 0.2, 0.0 ± 0.2],
        σk=[2.611 ± 0.076, 2.52 ± 0.077, 2.52 ± 0.077, 2.55 ± 0.092]
    )
end
