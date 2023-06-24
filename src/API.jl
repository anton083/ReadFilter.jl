
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

    writer = FASTAWriter(open(output_path, "w"))
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
        hits_d = max_scores_indices_d[findall(s -> s > score_threshold, max_scores_d)]
        hits_d = filter(idx -> idx[2] <= num_new_reads, hits_d)

        subref_indices_d = getindex.(hits_d, 1)
        read_indices_d = getindex.(hits_d, 2)
        hits_scores_d = vec(scores_d[hits_d])
        hits_byte_matrix_d = bases_to_bytes((reads_base_matrix_d[read_indices_d, :]))

        global_read_indices_d = read_indices_d .+ global_index_offset

        write_matched_reads(
            writer, hits_byte_matrix_d, global_read_indices_d, subref_indices_d, hits_scores_d)

        append!(flagged_reads, global_read_indices)

        n = length(flagged_reads)
        println("$n/$read_count ($(round(100*n/read_count, digits=2))%)")
    end
    close(reader)
    close(writer)

    #filter!(idx -> (idx <= read_count), flagged_reads)

    flagged_reads
end

    
# TODO: ask Kenta to make some stream that streams reads directly to byte matrix
# TODO: wrapper for loading a serialized reference matrix
# TODO: add check for homopolymers