
"""
Takes a sequence lengths and splits it into `subranges` of length `sublength` that overlap by `overlap`.
The overlap may be larger between the last and next to last subrange in case
"""
function get_subranges(full_length::Int, sublength::Int, overlap::Int)
    if full_length < sublength
        throw(ErrorException(
            "`full_length` ($full_length) needs to be greater than or equal to `sublength` ($sublength)"))
    end

    full_range = 1:full_length
    start_points = full_range.start:(sublength-overlap):(full_range.stop-L+1)
    subranges = [start:(start+sublength-1) for start in start_points]

    if subranges[end].stop < full_length
        push!(subranges, (full_length-sublength+1):full_length)
    end
    subranges
end

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
