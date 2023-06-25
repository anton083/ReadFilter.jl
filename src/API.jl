
function find_reads_gpu(
    ref_path::String,
    dataset_path::String,
    pident::Float64 = 0.9;
    k::Int = 6,
    subref_length::Int = 1024,
    read_chunk_size::Int = 100000,
    output_path::String = "filtered.fasta",
    check_alignments::Bool = false,
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
        n = length(flagged_reads)
        println("$n/$read_count ($(round(100*n/read_count, digits=2))%)")

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

        max_scores_indices_d = vec(CUDA.argmax(scores_d, dims=1))
        max_scores_d = scores_d[max_scores_indices_d]
        read_indices_d = findall(s -> s > score_threshold, max_scores_d)
        if !isempty(read_indices_d)
            read_indices_d = filter(idx -> idx <= num_new_reads, read_indices_d)
        end

        isempty(read_indices_d) && continue

        # Cartesian coordinates for score matrix
        match_cart_inds_d = max_scores_indices_d[read_indices_d]

        matched_reads_byte_matrix = Matrix{UInt8}(bases_to_bytes((reads_base_matrix_d[read_indices_d, :])))

        global_read_indices = Vector(read_indices_d .+ global_index_offset)
        subref_indices = Vector(getindex.(match_cart_inds_d, 1))
        matched_subrefs = subrefs[subref_indices]
        match_scores = Vector(max_scores_d[read_indices_d])

        reads = recreate_reads(matched_reads_byte_matrix, global_read_indices)
        read_matches = get_matches(reads, matched_subrefs, match_scores)

        if check_alignments
            assign_alignment_scores(read_matches)
            println(mean([rm.alignment_score for rm in read_matches]))
            filter!(rm -> rm.alignment_score > read_length / 2, read_matches)
            println(mean([rm.alignment_score for rm in read_matches]))
        end

        write_matches(writer, read_matches)

        append!(flagged_reads, [rm.read.idx for rm in read_matches])
    end
    close(reader)
    close(writer)

    flagged_reads
end


# TODO: ask Kenta to make some stream that streams reads directly to byte matrix
# TODO: wrapper for loading a serialized reference matrix
# TODO: add check for homopolymers
# TODO: revcomp! if matched to second half of subref vector (revcomps)
# TODO: don't count reads until after alignment filtering?