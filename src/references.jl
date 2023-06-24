
const Subrange = UnitRange{<:Integer}

struct Reference
    description::AbstractString
    sequence::LongDNA{4}
    length::Integer
    index::Integer
end

@inline Base.length(ref::Reference) = ref.length

function references(path::String, num_refs::Union{Integer, Float64} = Inf)
    records = read_records(path, num_refs)
    descs = description.(records)
    seqs = degap.(sequence.(LongDNA{4}, records))
    [Reference(desc, seq, length(seq), index) for (index, (desc, seq)) in enumerate(zip(descs, seqs))]
end


struct Subreference
    reference::Reference
    subrange::Subrange
    revcomp::Bool
end

@inline Base.length(subref::Subreference) = length(subref.subrange)

function subreferences(ref::Reference, subranges::Vector{<:Subrange})
    [Subreference(ref, subrange, false) for subrange in subranges]
end

revcomp(subref::Subreference) = Subreference(subref.reference, subref.subrange, !subref.revcomp)

function subreferences(
    refs::Vector{Reference},
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

    subranges_vec = get_subranges.(ref_lengths, sublength, read_length)
    subrefs = reduce(vcat, [subreferences(ref, subranges) for (ref, subranges) in zip(refs, subranges_vec)])

    # use divrem(i, N) to get (is_revcomp, j)
    append!(subrefs, revcomp.(subrefs))

    subrefs
end

function subreferences(
    ref_path::String,
    sublength::Integer,
    read_length::Integer,
)
    subreferences(references(ref_path), sublength, read_length)
end

function get_sequence(subref::Subreference)
    subseq = subref.reference.sequence[subref.subrange]
    subref.revcomp ? reverse_complement!(subseq) : subseq
end


function subref_kmer_matrix(
    subrefs::Vector{Subreference},
    subref_length::Integer,
    k::Integer,
)
    num_subrefs = length(subrefs)

    subref_byte_matrix_h = byte_matrix(num_subrefs, subref_length)
    for (j, subref) in enumerate(subrefs)
        byte_seq = codeunits(String(get_sequence(subref)))
        byte_seq_to_byte_matrix!(subref_byte_matrix_h, byte_seq, subref_length, j)
    end
    subref_base_matrix_d = subref_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases

    subref_kmer_matrix_d = kmer_count.GPU.row_bins(num_subrefs, k)
    kmer_count.GPU.kmer_count_rows!(subref_kmer_matrix_d, subref_base_matrix_d)

    subref_kmer_matrix_d
end