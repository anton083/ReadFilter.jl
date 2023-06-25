
struct Read
    seq::LongDNA{4}
    idx::Int
end

@inline Base.length(read::Read) = length(read.seq)

function recreate_reads(
    reads_byte_matrix::Matrix{UInt8},
    indices::Vector{Int},
)
    n = length(indices)
    reads = Vector{Match}(undef, n)
    rows = eachrow(reads_byte_matrix)
    for (i, (row, idx)) in enumerate(zip(rows, indices))
        seq = LongDNA{4}(String(row))
        read = Read(seq, idx)
        reads[i] = read
    end
    reads
end

struct Match
    read::Read
    subref::Subreference
    kmer_count_score::BinType
    alignment_score::Union{Missing, AlignScoreType}
end

function get_matches(
    reads::Vector{Read},
    subrefs::Vector{Subreference},
    kmer_count_scores::Vector{BinType},
)
    n = length(reads)
    matches = Vector{Match}(undef, n)
    for (i, read, subref, score) in enumerate(zip(reads, subrefs, kmer_count_scores))
        matches[i] = Match(read, subref, score, missing)
    end
    matches
end

function assign_alignment_scores(matches::Vector{Match})
    match = matches[1]
    params = Params(length(match.read), length(match.subref), 1, 1)

    for match in matches
        match.alignment_score = SWG_score(params, match.read.seq, get_sequence(match.subref))
    end
end