### ReadFilter.jl #

module ReadFilter

    using BioSequences
    using FASTX
    using LinearAlgebra
    using StatsBase
    using ProgressBars
    using StructArrays

    include("mutation.jl")
    include("kmers.jl")
    include("score.jl")
    include("io.jl")
    include("filter.jl")

    export
        kmer_count,
        kmer_count!,
        mutate,
        mutate!,
        score,
        estimate_score_threshold,
        read_records,
        single_filter_fasta

end
