### ReadFilter.jl #

module ReadFilter

    using BioSequences
    using FASTX
    using LinearAlgebra
    using StatsBase
    using StructArrays
    using CUDA

    include("utils.jl")
    include("io.jl")
    include("kmer_count/kmer_count.jl")
    include("references.jl")
    #include("score.jl")
    #include("mutation.jl")
    include("filter.jl")

    export filter_fasta_gpu

end