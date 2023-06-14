
module CPU

@inline kmer_count(k::Integer) = zeros(UInt16, 1 << 2k)

function kmer_count(seq::LongDNA{2}, k::Integer)
    bins = kmer_count(k)
    len = length(seq)
    mask = UInt(1 << 2k - 1)
    kmer = UInt(0)
    i = 0
    for data_int in seq.data
        for j in 0:2:62
            i += 1
            i > len && break
            kmer = kmer << 2 & mask + data_int >> j & 0b11
            @inbounds bins[kmer + 1] = k <= i
        end
    end
    bins
end

function kmer_count(seq::LongDNA{4}, k::Integer)
    bins = kmer_count(k)
    len = length(seq)
    mask = UInt(1 << 2k - 1)
    kmer = UInt(0)
    i = 0
    for data_int in seq.data
        for j in 0:4:60
            i += 1
            i > len && break
            kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
            bins[kmer + 1] = k <= i
        end
    end
    bins
end

function kmer_count!(bins::Vector{UInt16}, seq::LongDNA{2}, k::Integer)
    len = length(seq)
    mask = UInt(1 << 2k - 1)
    kmer = UInt(0)
    i = 0
    for data_int in seq.data
        for j in 0:2:62
            i += 1
            i > len && break
            kmer = kmer << 2 & mask + data_int >> j & 0b11
            @inbounds bins[kmer + 1] = k <= i
        end
    end
    bins
end

function kmer_count!(bins::Vector{UInt16}, seq::LongDNA{4}, k::Integer)
    len = length(seq)
    mask = UInt(1 << 2k - 1)
    kmer = UInt(0)
    i = 0
    for data_int in seq.data
        for j in 0:4:60
            i += 1
            i > len && break
            kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
            @inbounds bins[kmer + 1] = k <= i
        end
    end
    bins
end


function reference_kmer_matrix(refs::Vector{LongDNA{4}}, k::Integer)
    N = lastindex(refs)
    bin_matrix = zeros(UInt16, (N, 4^k))
    mask = UInt(4^k - 1)
    for (row, seq) in enumerate(refs)
        len = length(seq)
        kmer = UInt(0)
        i = 0
        for data_int in seq.data
            for j in 0:4:60
                i += 1
                i > len && break
                kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
                @inbounds bin_matrix[row,kmer+1] = k <= i
            end
        end
    end
    bin_matrix
end

function reference_kmer_matrix!(bin_matrix::Matrix{UInt16}, refs::Vector{LongDNA{4}}, k::Integer)
    fill!(bin_matrix, zero(eltype(bin_matrix)))
    mask = UInt(4^k - 1)
    for (row, seq) in enumerate(refs)
        len = length(seq)
        kmer = UInt(0)
        i = 0
        for data_int in seq.data
            for j in 0:4:60
                i += 1
                i > len && break
                kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
                @inbounds bin_matrix[row,kmer+1] = k <= i
            end
        end
    end
    bin_matrix
end


function read_kmer_matrix(reads::Vector{LongDNA{4}}, k::Integer)
    N = lastindex(reads)
    bin_matrix = zeros(UInt16, (4^k, N))
    mask = UInt(4^k - 1)
    for (column, seq) in enumerate(reads)
        len = length(seq)
        kmer = UInt(0)
        i = 0
        for data_int in seq.data
            for j in 0:4:60
                i += 1
                i > len && break
                kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
                @inbounds bin_matrix[kmer+1,column] = k <= i
            end
        end
    end
    bin_matrix
end

function read_kmer_matrix!(bin_matrix::Matrix{UInt16}, reads::Vector{LongDNA{4}}, k::Integer)
    fill!(bin_matrix, zero(eltype(bin_matrix)))
    mask = UInt(4^k - 1)
    for (column, seq) in enumerate(reads)
        len = length(seq)
        kmer = UInt(0)
        i = 0
        for data_int in seq.data
            for j in 0:4:60
                i += 1
                i > len && break
                kmer = kmer << 2 & mask + trailing_zeros(data_int >> j & 0b1111)
                @inbounds bin_matrix[kmer+1,column] = k <= i
            end
        end
    end
    bin_matrix
end

end