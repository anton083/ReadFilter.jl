
module kmer_count

    const BinType = Float16

    include("CPU.jl")
    include("GPU.jl")

    export CPU, GPU

end