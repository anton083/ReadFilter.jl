
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

function filter_fasta(
    read_index_set::Set{Int},
    fasta_path::String,
    output_file::String,
)
    reader = FASTAReader(open.(fasta_path, "r"))
    writer = FASTAWriter(open(output_file, "w"))

    for (read_index, record) in enumerate(reader)
        if read_index in read_index_set
            write(writer, record)
        end
    end

    close(reader)
    close(writer)
end

byte_matrix(num_seqs::Integer, seq_length::Integer) = zeros(UInt8, (num_seqs, seq_length))

"""
Does the bare minimum that needs to be done to the sequences before moving them to GPU.
"""



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

const BYTE_VECTOR = CuVector{UInt8}(UInt8.(collect("ACTG")))

@inline bases_to_bytes(base_matrix::CuMatrix{UInt8}) = BYTE_VECTOR[base_matrix .+ 0x01]


function longest_read_fasta(
    fasta_path::String,
    num_records::Int = 10,
)
    maximum(seqsize.(read_records(fasta_path, num_records)))
end

function get_indices_of_matches(
    scores_d::CuMatrix{BinType},
    score_thresholds_d::CuVector{BinType},
)
    match_bools_d = CUDA.reduce(max, scores_d .- score_thresholds_d .> 0, dims=1) # [1:(read_count - 1) % read_chunk_size + 1]
    indices_of_matches = findall(Vector(vec(match_bools_d)))
    indices_of_matches
end

# TODO: create a template LongDNA{2} of length `read_length` and fill data field with 2-bit bases?
function write_matched_reads(
    output_path::String,
    seq_matrix::Matrix{UInt8},
    read_indices::Vector{<:Integer},
    subref_indices::Vector{<:Integer},
    match_scores::Vector{BinType},
)
    writer = FASTAWriter(open(output_path, "w"))
    for (i, (read_idx, subref_idx, score)) in enumerate(zip(read_indices, subref_indices, match_scores))
        desc = "$read_idx $(subref_idx) $(round(score, digits=1))"
        seq = String(seq_matrix[i, :])
        write(writer, FASTARecord(desc, seq))
    end
    close(writer)
end