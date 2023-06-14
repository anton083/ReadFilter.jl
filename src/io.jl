
function read_records(fasta_path::String, num_records::Union{Integer, Float64} = Inf)
    num_records > 100000 && @warn "We might be loading a shit ton of records into memory"
    records = FASTARecord[]
    FASTAReader(open(fasta_path), copy=true) do reader
        i = 0
        for record in reader
            i += 1
            i > num_records && break
            push!(records, record)
        end
        i < num_records && !isinf(num_records) && @warn "Could only read $i records from `$fasta_path`"
    end
    records
end

byte_matrix(num_seqs::Integer, seq_length::Integer) = zeros(UInt8, (num_seqs, seq_length))

"""
Does the bare minimum that needs to be done to the sequences before moving them to GPU.
"""
function strings_to_byte_matrix(seqs::Vector{String})
    lengths = lastindex.(strings)
    @assert length(Set(lengths)) == 1 "Strings have different lengths"
    
    byte_vectors = Vector{UInt8}.(strings)
    byte_matrix = reduce(hcat, byte_vectors)
    byte_matrix
end

function strings_to_byte_matrix!(byte_matrix::Matrix{UInt8}, strings::Vector{String})
    lengths = lastindex.(strings) # lastindex is faster than length
    @assert length(Set(lengths)) == 1 "Strings have different lengths"
    @assert size(byte_matrix) == (length(strings), lengths[1]) "Sequences don't match byte_matrix. Byte matrix size $(size(byte_matrix)) doesn't match to $((length(seqs), lengths[1]))"
    
    for (i, str) in enumerate(strings)
        byte_matrix[i, :] = codeunits(str)
    end
    byte_matrix
end

function string_to_byte_matrix!(byte_matrix::Matrix{UInt8}, str::String, i::Integer)
    @assert size(byte_matrix, 2) == lastindex(str) "Sequence length doesn't match byte_matrix."
    byte_matrix[i, :] = codeunits(str)
end

#function load_seqs

"""
Uses bit manipulation to get convert ACGT into 4 unique bytes between 0 and 3.
Only designed for ACGT
For alphabetical order (both upper/lowercase) use: (byte - 0x01 - (byte % 0x20 == 0x03)) & 0x03
"""
@inline bytes_to_bases(byte_matrix::CuMatrix{UInt8}) = byte_matrix .>> 1 .& 0x03