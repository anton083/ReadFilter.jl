
module GPU

import ..kmer_count: BinType

using CUDA

# TODO: spaced_kmer_count (11011011 and 1001001001001001 instead of 111111)

export
    row_bins, kmer_count_rows!,
    column_bins, kmer_count_columns!

row_bins(N::Integer, k::Integer) = CUDA.zeros(BinType, (N, 4^k))

"""
The row version of bins has size `N x 4^k` 
The sequence matrix always has size `N x L`
"""
function kmer_count_rows!(bins::CuMatrix{BinType}, sequences::CuMatrix{UInt8}, k::Int)
    num_sequences, seq_len = size(sequences)
    CUDA.fill!(bins, zero(BinType))

    function kernel(sequences, bins, k, mask, num_sequences, seq_len)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if seq_idx <= num_sequences
            kmer = unsigned(0)
            for i in 1:k-1
                base = sequences[seq_idx, i]
                kmer = (kmer << 2) + base
            end
            for i in k:seq_len
                base = sequences[seq_idx, i]
                kmer = ((kmer << 2) & mask) + base
                CUDA.@allowscalar bins[seq_idx, kmer + 1] = one(BinType)
            end
        end
        return
    end

    mask = unsigned(4^k - 1)
    seq_len = size(sequences, 2)
    threads = 256
    blocks = ceil(Int, num_sequences / threads)

    @cuda threads=threads blocks=blocks kernel(sequences, bins, k, mask, num_sequences, seq_len)

    bins
end


column_bins(N::Integer, k::Integer) = CUDA.zeros(BinType, (4^k, N))

"""
The column version of bins has size `4^k x N` 
The sequence matrix always has size `N x L`
"""
function kmer_count_columns!(bins::CuMatrix{BinType}, sequences::CuMatrix{UInt8}, k::Int)
    num_sequences, seq_len = size(sequences)
    CUDA.fill!(bins, zero(BinType))

    function kernel(sequences, bins, k, mask, num_sequences, seq_len)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if seq_idx <= num_sequences
            kmer = unsigned(0)
            for i in 1:k-1
                base = sequences[seq_idx, i]
                kmer = (kmer << 2) + base
            end
            for i in k:seq_len
                base = sequences[seq_idx, i]
                kmer = ((kmer << 2) & mask) + base
                CUDA.@allowscalar bins[kmer + 1, seq_idx] = one(BinType)
            end
        end
        return
    end

    mask = unsigned(4^k - 1)
    seq_len = size(sequences, 2)
    threads = 256
    blocks = ceil(Int, num_sequences / threads)

    @cuda threads=threads blocks=blocks kernel(sequences, bins, k, mask, num_sequences, seq_len)

    bins
end

end