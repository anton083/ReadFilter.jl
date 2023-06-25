

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
        D[1, j] = Q[1, j] = -params.gap_open - (j - 2) * params.gap_extend
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

function SWG_score(seq1::LongDNA, seq2::LongDNA, gap_open::Int, gap_extend::Int)
    m, n = length(seq1), length(seq2)
    params = AlignParams(length(m), length(n), gap_open, gap_extend)
    SWG_score(params, seq1, seq2, gap_open, gap_extend)
end

function SWG_score(params::AlignParams, seq1::LongDNA, seq2::LongDNA, gap_open::Int, gap_extend::Int)
    m, n = length(seq1), length(seq2)
    @assert size(params.D) == (m+1, n+1)
    reset_matrices!(params)

    D, P, Q = params.D, params.P, params.Q
    
    for i in 2:(m+1)
        for j in 2:(n+1)
            P[i, j] = max(D[i-1, j] - gap_open, P[i-1, j] - gap_extend)
            Q[i, j] = max(D[i, j-1] - gap_open, Q[i, j-1] - gap_extend)
            match_score = seq1[i-1] == seq2[j-1] ? 1 : -1  # Assume match score of 1 and mismatch score of -1
            D[i, j] = max(0, D[i-1, j-1] + match_score, P[i, j], Q[i, j])
        end
    end
    
    maximum(D)
end
