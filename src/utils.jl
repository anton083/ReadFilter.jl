
degap(seq::LongDNA{4}) = filter(c -> c != DNA_Gap, seq)

max_in_columns(matrix::Array{Int, 2}) = mapslices(maximum, matrix, dims=1)

"""
Takes a sequence lengths and splits it into `subranges` of length `sublength` that overlap by `overlap`.
The overlap may be larger between the last and next to last subrange in case the last one doesn't perfectly fit.
"""
function get_subranges(full_length::Int, sublength::Int, overlap::Int)
    if full_length < sublength
        throw(ErrorException(
            "`full_length` ($full_length) needs to be greater than or equal to `sublength` ($sublength)"))
    end

    full_range = 1:full_length
    start_points = full_range.start:(sublength-overlap):(full_range.stop-sublength+1)
    subranges = [start:(start+sublength-1) for start in start_points]

    if subranges[end].stop < full_length
        push!(subranges, (full_length-sublength+1):full_length)
    end

    subranges
end

function random_subrange(L::Integer, R::Integer)
    start = rand(1:L-R+1)
    start:start+R-1
end
