### ReadFilter.jl #

module ReadFilter

    using BioSequences
    using FASTX
    using LinearAlgebra
    using StatsBase
    using CUDA
    using Plots

    include("utils.jl")
    include("io.jl")
    include("kmer_count/kmer_count.jl")
    include("references.jl")
    include("mutation.jl")
    include("simulation.jl")
    include("filter.jl")

    export filter_fasta_gpu

end