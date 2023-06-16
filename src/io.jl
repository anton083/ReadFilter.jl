
function read_records(fasta_path::String, num_records::Union{Integer, Float64} = Inf)
    100000 < num_records < Inf && @warn "We might be loading a shit ton of records into memory"
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
function strings_to_byte_matrix(strings::Vector{String})
    lengths = lastindex.(strings)
    @assert length(Set(lengths)) == 1 "Strings have different lengths"
    
    byte_vectors = Vector{UInt8}.(strings)
    byte_matrix = reduce(hcat, byte_vectors)
    byte_matrix
end


"""
The type returned by calling `codeunits(sequence(String, record))`.

It would be a little faster -- although less flexible -- to use `codeunits(sequence(record)).
"""
const CodeUnitsStringView = Base.CodeUnits{UInt8, String}

function byte_seq_to_byte_matrix!(
    byte_matrix::Matrix{UInt8},
    byte_seq::CodeUnitsStringView,
    len::Integer,
    i::Integer,
)
    new_len = min(len, size(byte_matrix, 2))
    byte_matrix[i, 1:new_len] = view(byte_seq, 1:new_len)
end

#function load_seqs

"""
Uses bit manipulation to get convert ACGT into 4 unique bytes between 0 and 3.
Only designed for ACGT
For alphabetical order (both upper/lowercase) use: (byte - 0x01 - (byte % 0x20 == 0x03)) & 0x03
"""
@inline bytes_to_bases(byte_matrix::CuMatrix{UInt8}) = byte_matrix .>> 1 .& 0x03