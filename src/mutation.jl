
const MUTATION_DICT = Dict{DNA, Tuple{DNA, DNA, DNA}}(
    DNA_A => (DNA_C, DNA_G, DNA_T),
    DNA_C => (DNA_A, DNA_G, DNA_T),
    DNA_G => (DNA_A, DNA_C, DNA_T),
    DNA_T => (DNA_A, DNA_C, DNA_G),
)

function mutate!(seq::LongDNA{4}, mut_rate::Float64) 
    for i in eachindex(seq)
        if rand() < mut_rate
            seq[i] = rand(get(MUTATION_DICT, seq[i], (DNA_A,)))
        end
    end
    seq
end

function mutate(seq::LongDNA{4}, mut_rate::Float64) 
    mutate!(deepcopy(seq), mut_rate)
end

function mutate!(seq::LongDNA{4}, mut_count::Int) 
    mutation_positions = sample(1:seq.len, mut_count, replace=false)
    println(mut_count)
    for mut_pos in mutation_positions
        seq[mut_pos] = rand(get(MUTATION_DICT, seq[mut_pos], (DNA_A,)))
    end
    seq
end

function mutate(seq::LongDNA{4}, mut_count::Int) 
    mutate!(deepcopy(seq), mut_count)
end