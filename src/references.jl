
# subsequences aren't views because we need the data field to be separate
@inline subsequence(seq::LongDNA{4}, subrange::UnitRange{Int}) = seq[subrange]

function subsequences(seq::LongDNA{4}, subranges::Vector{UnitRange{Int}})
    [subsequence(seq, subrange) for subrange in subranges]
end

function get_refs(ref_path::String, num_refs::Union{Integer, Float64} = Inf)
    degap.(sequence.(LongDNA{4}, read_records(ref_path, num_refs)))
end

"""
Takes a vector of reference sequences and splits them into equally-sized chunks (subreferences).
"""
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

    # use divrem(i, N) to get (is_revcomp, j)
    append!(subrefs, reverse_complement.(subrefs))

    subrefs
end


function subref_kmer_matrix(
    ref_path::String,
    subref_length::Integer,
    read_length::Integer,
    k::Integer,
)
    refs = get_refs(ref_path)
    subrefs = get_subrefs(refs, subref_length, read_length)
    num_subrefs = length(subrefs)

    subref_byte_matrix_h = byte_matrix(num_subrefs, subref_length)
    for (j, subref) in enumerate(subrefs)
        byte_seq = codeunits(String(subref))
        byte_seq_to_byte_matrix!(subref_byte_matrix_h, byte_seq, subref_length, j)
    end
    subref_base_matrix_d = subref_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases

    subref_kmer_matrix_d = kmer_count.GPU.row_bins(num_subrefs, k)
    kmer_count.GPU.kmer_count_rows!(subref_kmer_matrix_d, subref_base_matrix_d, k)

    subrefs, subref_kmer_matrix_d
end