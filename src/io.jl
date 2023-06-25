
function read_records(fasta_path::String, num_records::Union{Int, Float64} = Inf)
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

byte_matrix(num_seqs::Int, seq_length::Int) = zeros(UInt8, (num_seqs, seq_length))

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
    len::Int,
    i::Int,
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
bytes_to_bases(byte_matrix::CuMatrix{UInt8}) = byte_matrix .>> 1 .& 0x03

const CHAR_VECTOR = collect("ACTG")
const BYTE_VECTOR = UInt8.(CHAR_VECTOR)
const BYTE_VECTOR_d = cu(BYTE_VECTOR)

bases_to_bytes(base_matrix::Matrix{UInt8}) = BYTE_VECTOR[base_matrix .+ 0x01]
bases_to_bytes(base_matrix_d::CuMatrix{UInt8}) = BYTE_VECTOR_d[base_matrix_d .+ 0x01]

function longest_read_fasta(
    fasta_path::String,
    num_records::Int = 10,
)
    maximum(seqsize.(read_records(fasta_path, num_records)))
end

# TODO: create a template LongDNA{2} of length `read_length` and fill data field with 2-bit bases?
function write_matched_reads(
    writer::FASTAWriter,
    seq_matrix::Matrix{UInt8},
    read_indices::AbstractVector{Int},
    subref_indices::AbstractVector{Int},
    match_scores::AbstractVector{BinType},
)
    rows = eachrow(seq_matrix)
    for (seq, read_idx, subref_idx, score) in zip(rows, read_indices, subref_indices, match_scores)
        desc = "$read_idx $(subref_idx) $(round(score, digits=1))"
        seq = String(seq)
        write(writer, FASTARecord(desc, seq))
    end
end