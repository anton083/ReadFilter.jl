

const AlignScoreType = Float16
const MatrixType = Matrix{AlignScoreType}
const NegInf = AlignScoreType(-Inf)

struct AlignParams
    D::MatrixType
    P::MatrixType
    Q::MatrixType
    gap_open::Int
    gap_extend::Int
end

function reset_matrices!(params::AlignParams)
    fill!(params.D, NegInf)
    fill!(params.P, NegInf)
    fill!(params.Q, NegInf)

    params.D[1, 1] = 0
    for i in 2:size(params.D, 1)
        params.D[i, 1] = params.P[i, 1] = -params.gap_open - (i - 2) * params.gap_extend
    end
    for j in 2:size(params.D, 2)
        params.D[1, j] = params.Q[1, j] = -params.gap_open - (j - 2) * params.gap_extend
    end

    params
end

function AlignParams(m::Int, n::Int, gap_open::Int, gap_extend::Int)
    params = AlignParams(
        fill(NegInf, (m+1, n+1)),
        fill(NegInf, (m+1, n+1)),
        fill(NegInf, (m+1, n+1)),
        gap_open,
        gap_extend
    )
    reset_matrices!(params)
end

function SWG_align(seq1::LongDNA{4}, seq2::LongDNA{4}, gap_open::Int, gap_extend::Int)
    m, n = length(seq1), length(seq2)
    params = AlignParams(m, n, gap_open, gap_extend)
    SWG_align(seq1, seq2, params)
end

function SWG_align(seq1::LongDNA{4}, seq2::LongDNA{4}, params::AlignParams)
    m, n = length(seq1), length(seq2)
    @assert size(params.D) == (m+1, n+1)
    reset_matrices!(params)

    D, P, Q = params.D, params.P, params.Q
    gap_open, gap_extend = params.gap_open, params.gap_extend

    for i in 2:(m+1)
        for j in 2:(n+1)
            P[i, j] = max(D[i-1, j] - gap_open, P[i-1, j] - gap_extend)
            Q[i, j] = max(D[i, j-1] - gap_open, Q[i, j-1] - gap_extend)
            match_score = seq1[i-1] == seq2[j-1] ? 1 : -1  # Assume match score of 1 and mismatch score of -1
            D[i, j] = max(0, D[i-1, j-1] + match_score, P[i, j], Q[i, j])
        end
    end

    # Traceback from the maximum score in D
    alignment1 = []
    alignment2 = []
    I = argmax(D)
    i, j = I[1], I[2]

    while D[i, j] != 0
        if D[i, j] == D[i-1, j-1] + (seq1[i-1] == seq2[j-1] ? 1 : -1)
            push!(alignment1, seq1[i-1])
            push!(alignment2, seq2[j-1])
            i -= 1
            j -= 1
        elseif D[i, j] == P[i, j]
            push!(alignment1, seq1[i-1])
            push!(alignment2, '-')
            i -= 1
        else
            push!(alignment1, '-')
            push!(alignment2, seq2[j-1])
            j -= 1
        end
    end

    # Reverse the alignments as they are currently backwards
    alignment1 = reverse(alignment1)
    alignment2 = reverse(alignment2)

    return maximum(D), alignment1, alignment2
end