using TypedTables


function add_partition_info(events::Table, partitions::Table)::Table
    # NOTE: must preserve events ordering
    # FIXME: I'm sure the following is far from the best way to do it...
    rows = []
    for ev in events
        row = partitions[
            ev.timestamp .∈ partitions.span .&&
            ev.detector .== partitions.detector
        ]
        push!(rows, row)
    end

    return Table(merge(columns(events), columns(vcat(rows...))))
end

function add_event_idxs(partitions::Table, events::Table)::Table
    indices = fill(Vector{Int32}(), length(partitions))

    for (i, part) in enumerate(partitions)
        sel = events.timestamp .∈ part.span .&& events.detector .== part.detector
        indices[i] = findall(isequal(true), sel)
    end

    return Table((columns(partitions)..., event_idxs=indices))
end
