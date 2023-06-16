
"""
Assigns a score threshold to each subreference based on the average scores of uniformly mutated simulated reads.
"""
function get_score_thresholds(
    subrefs::Vector{LongDNA{4}},
    subref_kmer_matrix_d::CuMatrix{kmer_count.BinType},
    pident_threshold::Float64,
    k::Integer,
    read_length::Int,
    samples_per_subref::Integer = 100,
)
    subref_kmer_matrix_h = Matrix(subref_kmer_matrix_d)
    subref_len = length(subrefs[1])
    mut_count = trunc(Int, (1 - pident_threshold) * read_length)

    read_kmer_count_bins_d = kmer_count.GPU.column_bins(samples_per_subref, k)
    reads_byte_matrix_h = byte_matrix(samples_per_subref, read_length)

    score_thresholds = zeros(kmer_count.BinType, length(subrefs))
    for (i, subref) in enumerate(subrefs)
        subref_kmer_count = subref_kmer_matrix_h[i:i, :]

        for j in 1:samples_per_subref
            read_range = random_subrange(subref_len, read_length)
            read = mutate!(subref[read_range], mut_count)
            byte_seq = codeunits(String(read))
            byte_seq_to_byte_matrix!(reads_byte_matrix_h, byte_seq, read_length, j)
        end
        reads_base_matrix_d = reads_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases
        kmer_count.GPU.kmer_count_columns!(reads_kmer_count_bins_d, reads_base_matrix_d, k)
        score_thresholds[i] = mean(CuMatrix(subref_kmer_count) * read_kmer_count_bins_d)
    end
    
    score_thresholds
end
