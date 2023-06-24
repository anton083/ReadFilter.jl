
"""
Assigns a score threshold to each subreference based on the average scores of uniformly mutated simulated reads.
"""
function get_score_thresholds(
    subrefs::Vector{Subreference},
    subref_kmer_matrix_d::CuMatrix{BinType},
    pident_threshold::Float64,
    k::Integer,
    subref_length::Integer,
    read_length::Integer,
    samples_per_subref::Integer = 100,
)
    mut_count = trunc(Int, (1 - pident_threshold) * read_length)

    reads_kmer_matrix_d = kmer_count.GPU.column_bins(samples_per_subref, k)
    reads_byte_matrix_h = byte_matrix(samples_per_subref, read_length)

    score_thresholds = zeros(BinType, length(subrefs))
    for (i, subref) in enumerate(subrefs)
        subref_kmer_count = subref_kmer_matrix_d[i:i, :]

        for j in 1:samples_per_subref
            read_range = random_subrange(subref_length, read_length)
            read = mutate!(get_sequence(subref)[read_range], mut_count)
            byte_seq = codeunits(String(read))
            byte_seq_to_byte_matrix!(reads_byte_matrix_h, byte_seq, read_length, j)
        end
        reads_base_matrix_d = bytes_to_bases(CuMatrix{UInt8}(reads_byte_matrix_h))
        kmer_count.GPU.kmer_count_columns!(reads_kmer_matrix_d, reads_base_matrix_d, k)
        score_thresholds[i] = BinType(mean(Float32, subref_kmer_count * reads_kmer_matrix_d))
    end
    
    score_thresholds
end
