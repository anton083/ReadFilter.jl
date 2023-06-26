
struct Read
    seq::LongDNA{4}
    idx::Int
end

@inline Base.length(read::Read) = length(read.seq)

function recreate_reads(
    reads_byte_matrix::Matrix{UInt8},
    indices::Vector{Int},
    global_indices::Vector{Int},
)
    n = length(indices)
    reads = Vector{Read}(undef, n)
    for (i, (idx, global_idx)) in enumerate(zip(indices, global_indices))
        seq = LongDNA{4}(reads_byte_matrix[idx, :])
        read = Read(seq, global_idx)
        reads[i] = read
    end
    reads
end

mutable struct Match
    read::Read
    subref::Subreference
    kmer_count_score::BinType
    alignment_score::Union{Missing, AlignScoreType}
end

function get_matches(
    reads::Vector{Read},
    subrefs::Vector{Subreference},
    all_subrefs::Vector{Subreference},
    kmer_count_scores::Vector{BinType},
)
    n = length(reads)
    matches = Vector{Match}(undef, n)
    for (i, (read, subref, score)) in enumerate(zip(reads, subrefs, kmer_count_scores))
        match = Match(read, subref, score, missing)
        matches[i] = match
        alignment_score, a1, a2 = SWG_align(read.seq, get_sequence(subref), 1, 1)
        if alignment_score > 40
            println(score)
            println(alignment_score)
            print_alignment(a1, a2)
            print_alignment(a2, a1)
        end
    end
    matches
end

function assign_alignment_scores(matches::Vector{Match})
    match = matches[1]
    params = AlignParams(length(match.read), length(match.subref), 1, 1)
    for match in matches
        match.alignment_score, a1, a2 = SWG_align(match.read.seq, get_sequence(match.subref), params)
        println(match.kmer_count_score)
        println(match.alignment_score)
        print_alignment(a1, a2)
        print_alignment(a2, a1)
    end
end

# TODO: keep track of subref index. match.subref_idx field would be easiest to implement
@inline function write_match(writer::FASTAWriter, match::Match)
    read_idx = match.read.idx
    ref_index = match.subref.reference.idx
    score1 = round(match.kmer_count_score, digits=1)
    score2 = round(match.alignment_score, digits=1)
    desc = "$read_idx r$(ref_index) $(score1):$(score2)"
    write(writer, FASTARecord(desc, match.read.seq))
end

# TODO: create a template LongDNA{2} of length `read_length` and fill data field with 2-bit bases?
function write_matches(writer::FASTAWriter, matches::Vector{Match})
    for match in matches
        write_match(writer, match)
    end
end