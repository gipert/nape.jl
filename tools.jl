using LegendHDF5IO
import IntervalSets
using TypedTables

function _make_events_columns()
    return (
        Vector{Int64}(),
        Vector{Symbol}(),
        Vector{Float32}()
    )
end

function _make_partitions_columns()
    return (
        Vector{IntervalSets.ClosedInterval{Int64}}(),
        Vector{Symbol}(),
        Vector{Float32}(),
        Vector{Measurement{Float32}}(),
        Vector{Measurement{Float32}}(),
        Vector{Measurement{Float32}}()
    )
end

# FIXME: I want this method to modify events (!), but how do I replace the
# "events" object? is it even possible?
function _add_partition_info(events::Table, partitions::Table)::Table
    # NOTE: must preserve events ordering
    # FIXME: I'm sure the following is far from the best way to do it...
    rows = []
    part_idx = Vector{Int32}()

    for ev in events
        sel = ev.timestamp .∈ partitions.span .&& ev.detector .== partitions.detector
        row = partitions[sel]
        length(row) != 1 && @error "Could not find partition event $ev belongs to or multiple partitions found"

        push!(part_idx, first(findall(isequal(true), sel)))
        push!(rows, row)
    end

    return Table(merge(columns(events), columns(vcat(rows...)), (part_idx=part_idx,)))
end

function _add_parameter_index(events::Table, partitions::Table)::Table
    return Table((
        columns(events)...,
        epar_idx=indexin(events.part_idx, unique(events.part_idx))
    ))
end

function _add_event_idxs(partitions::Table, events::Table)::Table
    indices = fill(Vector{Int32}(), length(partitions))

    for (i, part) in enumerate(partitions)
        sel = events.timestamp .∈ [part.span] .&& events.detector .== part.detector
        indices[i] = findall(isequal(true), sel)
    end

    return Table((columns(partitions)..., event_idxs=indices))
end

function get_data(experiment::Symbol, from_cache::Bool=false)::NamedTuple
    @info "getting data for $experiment..."
    if from_cache
        #...
        # return (events=events, partitions=partitions)
    elseif startswith(string(experiment), "gerdaI_")
        events = read_events_gerdaI(experiment)
        partitions = read_partitions_gerdaI()
    elseif experiment == :gerdaII
        events = read_events_gerdaII()
        partitions = read_partitions_gerdaII()
    elseif experiment == :legend200
        events = read_events_legend200()
        partitions = read_partitions_legend200()
    elseif startswith(string(experiment), "majorana_")
        events = read_events_majorana(experiment)
        partitions = read_partitions_majorana()
    else
        @error "experiment $experiment not known"
    end

    events = _add_partition_info(events, partitions)
    events = _add_parameter_index(events, partitions)
    partitions = _add_event_idxs(partitions, events)

    return (events=events, partitions=partitions)
end

function cache_data(events, partitions, experiment)
    store = lh5open(joinpath(@__DIR__), ".cache", "$(string(experiment)).lh5")
    # store["events"] = ...
end
