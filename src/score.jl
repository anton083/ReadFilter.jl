
function score(
    reference_kmer_count::Vector{UInt32},
    record_kmer_count::Vector{UInt32},
    read_length::Int,
)
    (reference_kmer_count ⋅ record_kmer_count) / read_length
end

function estimate_score_threshold(
    reference::LongDNA{4},
    pident_threshold::Float64;
    k = 7,
    sample_count = 1000,
)
    mut_count = floor(Int, (1 - pident_threshold) * length(reference))
    scores = Float64[]
    reference_kmer_count = kmer_count(reference, k)
    read_kmer_count = kmer_count(k)
    for i in 1:sample_count
        read = mutate(reference, mut_count)
        fill!(read_kmer_count, zero(UInt32))
        kmer_count!(read_kmer_count, read, k)
        push!(scores, score(reference_kmer_count, read_kmer_count, length(read)))
    end
    sum(scores) / sample_count, reference_kmer_count
end

# TODO: account for paired read layout (make read concat with revcomp?) 