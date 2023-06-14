
function score(
    reference_kmer_count::Vector{UInt16},
    record_kmer_count::Vector{UInt16},
    read_length::Int,
)
    (reference_kmer_count ⋅ record_kmer_count) / read_length
end

function estimate_score_threshold(
    reference::LongDNA{4},
    pident_threshold::Float64;
    k = 7,
    num_samples = 1000,
)
    mut_count = floor(Int, (1 - pident_threshold) * length(reference))
    scores = Float64[]
    reference_kmer_count = kmer_count(reference, k)
    read_kmer_count = kmer_count(k)
    for _ in 1:num_samples
        read = mutate(reference, mut_count)
        fill!(read_kmer_count, zero(UInt16))
        kmer_count!(read_kmer_count, read, k)
        push!(scores, score(reference_kmer_count, read_kmer_count, length(read)))
    end
    sum(scores) / num_samples, reference_kmer_count
end

function estimate_score_threshold2(
    reference::LongDNA{4},
    pident_threshold::Float64,
    read_length::Int = 100;
    k = 7,
    num_samples = 1000,
)
    ref_len = length(reference)
    mut_count = floor(Int, (1 - pident_threshold) * read_length)
    scores = Float64[]
    reference_kmer_count = kmer_count(reference, k)
    read_kmer_count = kmer_count(k)
    for _ in 1:num_samples
        read_start = rand(1:ref_len-read_length+1)
        read = mutate(reference[read_start:read_start+read_length-1], mut_count)
        fill!(read_kmer_count, zero(UInt16))
        kmer_count!(read_kmer_count, read, k)
        push!(scores, score(reference_kmer_count, read_kmer_count, read_length))
    end
    sum(scores) / num_samples, reference_kmer_count
end
