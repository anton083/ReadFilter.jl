
function find_reads_gpu(
    ref_path::String,
    dataset_path::String,
    pident::Float64 = 0.9;
    k::Integer = 6,
    subref_length::Integer = 1024,
    read_chunk_size::Integer = 100000,
    output_path::String = "filtered.fasta"
)
    read_length = longest_read_fasta(dataset_path)
    subrefs = subreferences(ref_path, subref_length, read_length)
    subref_length = length(subrefs[1])

    subref_kmer_matrix_d = subref_kmer_matrix(subrefs, subref_length, k)

    score_thresholds_d = CuVector(get_score_thresholds(
        subrefs, subref_kmer_matrix_d, pident, k, subref_length, read_length))

    score_threshold = BinType(mean(Float32, score_thresholds_d))    

    reads_kmer_matrix_d = kmer_count.GPU.column_bins(read_chunk_size, k)
    reads_byte_matrix_h = byte_matrix(read_chunk_size, read_length)

    # Pre-allocate scores_d? CUDA.zeros(BinType, (ref_count, read_chunk_size))
    flagged_reads = Int64[]

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
        num_new_reads = (read_count - 1) % read_chunk_size + 1
        global_index_offset = read_count - num_new_reads

        reads_base_matrix_d = bytes_to_bases(CuMatrix{UInt8}(reads_byte_matrix_h))
        kmer_count.GPU.kmer_count_columns!(reads_kmer_matrix_d, reads_base_matrix_d, k)
        scores_d = subref_kmer_matrix_d * reads_kmer_matrix_d

        max_scores_indices_d = CUDA.argmax(scores_d, dims=1)
        max_scores_d = scores_d[max_scores_indices_d]
        hits = Vector(max_scores_indices_d[findall(s -> s > score_threshold, max_scores_d)])
        filter!(idx -> idx[2] <= num_new_reads, hits)

        subref_indices = getindex.(hits, 1)
        read_indices = getindex.(hits, 2)
        hits_scores = Vector(vec(scores_d[hits]))
        hits_seqs = bases_to_bytes(Matrix{UInt8}(reads_base_matrix_d[read_indices, :]))

        write_matched_reads(
            output_path, hits_seqs,
            read_indices, subref_indices, hits_scores,
        )

        global_read_indices = read_indices .+ global_index_offset

        append!(flagged_reads, global_read_indices)

        n = length(flagged_reads)
        println("$n/$read_count ($(round(100*n/read_count, digits=2))%)")
    end
    close(reader)

    #filter!(idx -> (idx <= read_count), flagged_reads)

    flagged_reads
end

    
# TODO: ask Kenta to make some stream that streams reads directly to byte matrix
# TODO: wrapper for loading a serialized reference matrix
# TODO: add check for homopolymers