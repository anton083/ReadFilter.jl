
module kmer_count

    import ..ReadFilter: BinType

    #include("CPU.jl")
    include("GPU.jl")

    export
        #CPU,
        GPU

end