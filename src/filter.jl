

function chunked_filter_fasta(
    reference_path::String,
    dataset_path::String;
    pident::Float64 = 0.7,
    k::Integer = 7,
    subref_length::Integer = 2000,
    read_chunk_size::Integer = 10000,
)

    ref_seqs = degap.(sequence.(LongDNA{4}, read_records(reference_path)))
    ref_subranges = sequence_subranges.(length.(ref_seqs), subref_length, 100)
    subrefs = reduce(vcat, subseqs.(ref_seqs, ref_subranges))

    score_threshold, _ = estimate_score_threshold2(subrefs[1], pident, 90, k = k, sample_count = 1000)

    reference_bin_matrix = reference_kmer_matrix(subrefs, k)
    
    read_bin_matrix = zeros(UInt16, (4^k, read_chunk_size))
    read_chunk = StructVector{Read}(undef, read_chunk_size)

    matrix_product = zeros(Int, (length(subrefs), read_chunk_size))

    reader = FASTAReader(open(dataset_path), copy=false)
    num_read_total = 0
    hits = 0
    while !eof(reader)
        for (i, record) in enumerate(reader)
            num_read_total += 1
            read_chunk.seq[i] = sequence(LongDNA{4}, record)
            read_chunk.len[i] = seqsize(record)
            read_chunk.idx[i] = i
            i == read_chunk_size && break
        end
        read_kmer_matrix!(read_bin_matrix, read_chunk.seq, k)
        mul!(matrix_product, reference_bin_matrix, read_bin_matrix)

        hits += count(score -> score > score_threshold, (max_in_columns(matrix_product) / mean(read_chunk.len)))
        #if num_read_total % 10000 == 0 println(num_read_total) end
        println("$hits/$num_read_total")
    end
    close(reader)
end

function filter_fasta_gpu(
    ref_path::String,
    dataset_path::String,
    pident::Float64 = 0.9;
    k::Integer = 6,
    subref_length::Integer = 1024,
    read_chunk_size::Integer = 100000,
    read_length::Integer = 90,
    num_refs::Union{Integer, Float64} = Inf,
)
    refs = get_refs(ref_path::String, num_refs)
    subrefs = get_subrefs(refs, subref_length, read_length)
    num_subrefs = length(subrefs)

    subrefs_base_matrix_d = strings_to_byte_matrix(String.(subrefs)) |> CuMatrix{UInt8} |> bytes_to_bases
    subrefs_kmer_count_bins_d = kmer_count.GPU.row_bins(num_subrefs, k)
    kmer_count.GPU.kmer_count_rows!(subrefs_kmer_count_bins_d, subrefs_base_matrix_d, k)

    score_threshold, _ = estimate_score_threshold2(subrefs[1], pident, read_length, k = k, num_samples = 1000)
    
    reads_kmer_count_bins_d = kmer_count.GPU.column_bins(read_chunk_size, k)
    reads_byte_matrix_h = byte_matrix(read_chunk_size, read_length)
    
    # TODO: ask Kenta to make some stream that streams reads directly to byte matrix

    reader = FASTAReader(open(dataset_path), copy=false)
    num_read_total = 0
    hits = 0
    while !eof(reader)
        for (i, record) in enumerate(reader)
            num_read_total += 1
            str = sequence(String, record)
            string_to_byte_matrix!(reads_byte_matrix_h, str, i)
            i == read_chunk_size && break
        end
        reads_base_matrix_d = reads_byte_matrix_h |> CuMatrix{UInt8} |> bytes_to_bases
        kmer_count.GPU.kmer_count_columns!(reads_kmer_count_bins_d, reads_base_matrix_d, k)
        mul!(matrix_product, reference_bin_matrix, read_bin_matrix)

        hits += count(score -> score > score_threshold, (max_in_columns(matrix_product) / mean(read_chunk.len)))
    end
    close(reader)
end