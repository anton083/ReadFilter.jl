
# subsequences aren't views because we need the data field to be separate
@inline subsequence(seq::LongDNA{4}, subrange::UnitRange{Int}) = seq[subrange]

function subsequences(seq::LongDNA{4}, subranges::Vector{UnitRange{Int}})
    [subsequence(seq, subrange) for subrange in subranges]
end

function get_refs(ref_path::String, num_refs::Union{Integer, Float64} = Inf)
    degap.(sequence.(LongDNA{4}, read_records(ref_path, num_refs)))
end

function get_subrefs(
    refs::Vector{LongDNA{4}},
    sublength::Integer,
    read_length::Integer,
)
    ref_lengths = length.(refs)

    smallest_length = minimum(ref_lengths)
    if smallest_length < sublength
        @warn """Reference sequence shorter than `sublength` found. Switching to smaller sublength, $smallest_length..."*
                 This may result in more subreferences and longer runtime."""
        sublength = smallest_length
    end

    subranges = get_subranges.(ref_lengths, sublength, read_length)
    subrefs = reduce(vcat, subsequences.(refs, subranges))

    subrefs
end