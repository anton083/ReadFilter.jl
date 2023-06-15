
function filter_fasta_gpu(
    ref_path::String,
    dataset_path::String,
    pident::Float64 = 0.9;
    k::Integer = 6,
    subref_length::Integer = 1024,
    read_chunk_size::Integer = 100000,
    read_length::Integer = 90,
    num_refs::Union{Integer, Float64} = Inf,
)
    refs = get_refs(ref_path::String, num_refs)
    subrefs = get_subrefs(refs, subref_length, read_length)
    num_subrefs = length(subrefs)

    subrefs_base_matrix_d = strings_to_byte_matrix(String.(subrefs)) |> CuMatrix{UInt8} |> bytes_to_bases
    subrefs_kmer_count_bins_d = kmer_count.GPU.row_bins(num_subrefs, k)
    kmer_count.GPU.kmer_count_rows!(subrefs_kmer_count_bins_d, subrefs_base_matrix_d, k)

    #score_threshold, _ = estimate_score_threshold2(subrefs[1], pident, read_length, k = k, num_samples = 1000)
    
    reads_kmer_count_bins_d = kmer_count.GPU.column_bins(read_chunk_size, k)
    reads_byte_matrix_h = byte_matrix(read_chunk_size, read_length)
    
    # TODO: ask Kenta to make some stream that streams reads directly to byte matrix

    reader = FASTAReader(open(dataset_path), copy=false)
    num_reads_total = 0

    #scores_d = CUDA.zeros(kmer_count.BinType, (ref_count, read_chunk_size))
    all_max_scores = kmer_count.BinType[]

    while !eof(reader)
        for (i, record) in enumerate(reader)
            num_reads_total += 1
            str = sequence(String, record)
            string_to_byte_matrix!(reads_byte_matrix_h, str, i)
            i == read_chunk_size && break
        end
        reads_base_matrix_d = reads_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases
        kmer_count.GPU.kmer_count_columns!(reads_kmer_count_bins_d, reads_base_matrix_d, k)
        scores_d = subrefs_kmer_count_bins_d * reads_kmer_count_bins_d
        max_scores = Array(CUDA.reduce(max, scores_d, dims=1))

        append!(all_max_scores, max_scores)
        println(num_reads_total)
    end
    close(reader)

    all_max_scores
    #histogram(all_max_scores, bins=50)
    #savefig("plot$(time()).png")
end