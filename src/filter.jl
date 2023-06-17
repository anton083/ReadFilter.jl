
function filter_fasta_gpu(
    ref_path::String,
    dataset_path::String,
    pident::Float64 = 0.9;
    k::Integer = 6,
    subref_length::Integer = 1024,
    read_chunk_size::Integer = 100000,
    read_length::Integer = 90,
)
    subrefs, subref_kmer_matrix_d = subref_kmer_matrix(ref_path, subref_length, read_length, k)

    score_thresholds_d = CuVector(get_score_thresholds(subrefs, subref_kmer_matrix_d, pident, k, read_length))
    println(mean(score_thresholds_d))
    
    reads_kmer_matrix_d = kmer_count.GPU.column_bins(read_chunk_size, k)
    reads_byte_matrix_h = byte_matrix(read_chunk_size, read_length)

    # Pre-allocate scores_d? CUDA.zeros(kmer_count.BinType, (ref_count, read_chunk_size))
    flagged_reads = Int64[]
    all_max_scores = kmer_count.BinType[]

    reader = FASTAReader(open(dataset_path), copy=false)
    read_count = 0
    while !eof(reader)
        for (i, record) in enumerate(reader)
            read_count += 1
            byte_seq = codeunits(sequence(String, record))
            len = seqsize(record)
            byte_seq_to_byte_matrix!(reads_byte_matrix_h, byte_seq, len, i)
            i == read_chunk_size && break
        end
        reads_base_matrix_d = reads_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases
        kmer_count.GPU.kmer_count_columns!(reads_kmer_matrix_d, reads_base_matrix_d, k)
        scores_d = subref_kmer_matrix_d * reads_kmer_matrix_d

        match_bools_d = CUDA.reduce(max, scores_d .- score_thresholds_d .> 0, dims=1)
        match_indices = findall(Vector(vec(match_bools_d)))
        append!(match_indices, match_indices)

        max_scores = Array(CUDA.reduce(max, scores_d, dims=1))
        append!(all_max_scores, view(max_scores, 1:(read_count - 1) % read_chunk_size + 1))
        println("$(maximum(max_scores)), $(mean(max_scores))")
        println("$(length(flagged_reads))/$read_count")
    end
    close(reader)

    flagged_reads, all_max_scores
end

    
# TODO: ask Kenta to make some stream that streams reads directly to byte matrix
# TODO: wrapper for loading a serialized reference matrix