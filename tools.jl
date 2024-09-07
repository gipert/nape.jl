using LegendHDF5IO
using TypedTables


# FIXME: I want this method to modify events (!), but how do I replace the
# "events" object? is it even possible?
function _add_partition_info(events::Table, partitions::Table)::Table
    # NOTE: must preserve events ordering
    # FIXME: I'm sure the following is far from the best way to do it...
    rows = []; part_idx = []
    for ev in events
        sel = ev.timestamp .∈ partitions.span .&& ev.detector .== partitions.detector
        row = partitions[sel]
        length(row) != 1 && @error "Could not find partition event $ev belongs to or multiple partitions found"

        push!(part_idx, first(findall(isequal(true), sel)))
        push!(rows, row)
    end

    return Table(merge(columns(events), columns(vcat(rows...)), (part_idx=part_idx,)))
end

function _add_event_idxs(partitions::Table, events::Table)::Table
    indices = fill(Vector{Int32}(), length(partitions))

    for (i, part) in enumerate(partitions)
        sel = events.timestamp .∈ [part.span] .&& events.detector .== part.detector
        indices[i] = findall(isequal(true), sel)
    end

    return Table((columns(partitions)..., event_idxs=indices))
end

function get_data(experiment::Symbol, from_cache=false)::NamedTuple
    if from_cache
        #...
        # return (events=events, partitions=partitions)
    elseif experiment == :gerdaII
        events = read_events_gerdaII()
        partitions = read_partitions_gerdaII()
    elseif experiment == :legend200
        events = read_events_legend200()
        partitions = read_partitions_legend200()
    elseif string(experiment)[1:8] == "majorana"
        events = read_events_majorana(experiment)
        partitions = read_partitions_majorana()
    else
        @error "experiment $experiment not known"
    end

    events = _add_partition_info(events, partitions)
    partitions = _add_event_idxs(partitions, events)

    return (events=events, partitions=partitions)
end

function cache_data(events, partitions, experiment)
    store = lh5open(joinpath(@__DIR__), ".cache", "$(string(experiment)).lh5")
    # store["events"] = ...
end
