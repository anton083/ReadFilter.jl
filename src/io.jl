
function read_records(fasta_path::String, record_count::Integer)
    records = FASTARecord[]
    FASTAReader(open(fasta_path), copy=true) do reader
        i = 0
        for record in reader
            i += 1
            i > record_count && break
            push!(records, record)
        end
    end
    records
end

function read_records(fasta_path::String)
    reader = FASTAReader(open(fasta_path))
    records = collect(reader)
    close(reader)
    records
end


function sequence_subranges(seq_length::Int, L::Int, R::Int)
    seq_range = 1:seq_length
    start_points = seq_range.start:(L-R):(seq_range.stop-L+1)
    subranges = [s:(s+L-1) for s in start_points]

    if subranges[end].stop < seq_length
        push!(subranges, (seq_length-L+1):seq_length)
    end
    subranges
end

export sequence_subranges

# subseqs aren't views cause we need the seq.data to be separate
# well, actually, I guess we could read the subseqs from seq.data
# performance doesn't really matter here since there aren't a lot of refs
# it's linear time and space complexity regardless
@inline subseq(seq::LongDNA{4}, subrange::UnitRange{Int}) = seq[subrange]

subseqs(seq::LongDNA{4}, subranges::Vector{UnitRange{Int}}) = [subseq(seq, subrange) for subrange in subranges]