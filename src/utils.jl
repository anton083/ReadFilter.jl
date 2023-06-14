
degap(seq::LongDNA{4}) = filter(c -> c != DNA_Gap, seq)

max_in_columns(matrix::Array{Int, 2}) = mapslices(maximum, matrix, dims=1)