


function single_filter_fasta(
    reference_path::String,
    dataset_path::String;
    pident::Float64 = 0.7,
    k::Integer = 7,
)
    reference = sequence(LongDNA{4}, first(read_records(reference_path, 1)))
    score_threshold, reference_kmer_count = estimate_score_threshold(reference, pident, k = k, sample_count = 1000)

    reader = FASTAReader(open(dataset_path), copy=false)
    
    marked_records = Set(Int[])
    
    record_kmer_count = kmer_count(k)
    record_index = 0
    for record in ProgressBar(reader)
        record_index += 1
        fill!(record_kmer_count, zero(UInt32))
        kmer_count!(record_kmer_count, sequence(LongDNA{4}, record), k)
        record_score = score(reference_kmer_count, record_kmer_count, seqsize(record))
        record_score > score_threshold && push!(marked_records, record_index)
    end

    close(reader)

    return length(marked_records), marked_records
end

struct Read
    seq::LongDNA{4}
    len::Int64
    idx::Int64
end

function chunked_filter_fasta(
    reference_path::String,
    dataset_path::String;
    pident::Float64 = 0.7,
    k::Integer = 7,
    subref_length::Integer = 2048,
    read_chunk_size::Integer = 1000,
)

    ref_seqs = sequence.(LongDNA{4}, read_records(reference_path))
    ref_subranges = sequence_subranges.(length.(ref_seqs), subref_length, 100)
    subrefs = reduce(vcat, subseqs.(ref_seqs, ref_subranges))

    reference_bin_matrix = reference_kmer_matrix(subrefs, 7)

    read_bin_matrix = zeros(UInt16, (1000, 4^k))

    read_chunk = StructVector{Read}(undef, read_chunk_size)

    reader = FASTAReader(open(dataset_path), copy=false)
    read_count_total = 0
    while !eof(reader)
        for (i, record) in enumerate(reader)
            read_count_total += 1
            read_chunk.seq[i] = sequence(LongDNA{4}, record)
            read_chunk.len[i] = seqsize(record)
            read_chunk.idx[i] = i
            i == 100 && break
        end
        if read_count_total % read_chunk_size == 0
            read_kmer_matrix!(read_bin_matrix, read_chunk.seq, k)
            reference_bin_matrix * read_bin_matrix
        end
    end

    close(reader)
end

export chunked_filter_fasta