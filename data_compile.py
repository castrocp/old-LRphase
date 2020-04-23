import bamnostic as bs
import os
import vcf
import time
import re
import collections

def cigar_alignment_padded(seq=None, cigar=None, start_pos = 0, qualities=None, base_qual_thresh=0, query=False):
    """
    Note:
        Modified function from BAMnostic to pad for deletions.
        Any clipping results in the removal of those bases. If an insertion is seen in
        the CIGAR, those bases are removed from the sequence. If a deletion is seen in
        the CIGAR, those bases are padded with a period ('-') symbol.
    Args:
        seq (str): string sequence of the aligned segment.
        cigar (str): the cigar string or `cigartuple` of the aligned segment.
        start_pos (int): the first aligned position of the read
        qualities (:py:obj:`array.array`): base quality array from read
    Yields:
        (:py:obj:`tuple` of :py:obj:`str` and :py:obj:`int`): nucleotide base and index position of that base relative to reference
    """ 
    cigar = bs.utils.check_cigar_arg(cigar)
    cigar_aligned = ''
    algn_seg = {}
    last_cigar_pos = 0
    for op, n_ops in cigar:
        op_id = op if type(op) is int else op[1]
        if op_id == 5:  # BAM_CHARD_CLIP: skip hard clip CIGAR ops
            pass
        elif op_id in {1, 4}:  # BAM_CINS or BAM_CSOFT_CLIP: remove from sequence
            if query and op_id == 1:
                seg_seq = seq[last_cigar_pos:last_cigar_pos + n_ops]
                if qualities is not None:
                    seg_qual = qualities[last_cigar_pos:last_cigar_pos + n_ops]
                for index, base in enumerate(seg_seq):
                    if qualities is not None:
                        if seg_qual[index] >= base_qual_thresh:
                            yield base, start_pos
                    else:
                        yield base, start_pos
            last_cigar_pos += n_ops
        elif op_id == 3:  # BAM_CREF_SKIP: intron or large gaps
            start_pos += n_ops
        elif op_id == 2:  # BAM_CDEL: pad for deletions
            for i in range(n_ops):
                yield '-', start_pos
                start_pos += 1
            start_pos += n_ops
        elif op_id in {0, 7, 8}:  # matches (uses both sequence match & mismatch)
            seg_seq = seq[last_cigar_pos:last_cigar_pos + n_ops]
            if qualities is not None:
                seg_qual = qualities[last_cigar_pos:last_cigar_pos + n_ops]
            for index, base in enumerate(seg_seq):
                if qualities is not None:
                    if seg_qual[index] >= base_qual_thresh:
                        yield base, start_pos
                else:
                    yield base, start_pos
                start_pos += 1
            last_cigar_pos += n_ops
        else:
            raise ValueError('Invalid CIGAR string: {}'.format(op))


def process_data(output_dir, reference_genome, vcf_file):  
    start_process_time= time.time()
    print("************* Processing Reads")

    # Read alignments from the BAM file
    dirname = os.path.join(output_dir, "alignment-output")
    bam_file = bs.AlignmentFile(dirname + "/alignment.sorted.bam", 'rb') 

    # Data for all SNPs overlapping the aligned read is stored here and passed to data_analysis.py
    read_snps = []
    
# Count the number of reads being processed and filtered
    total_alignment_count = 0
    query_alignment_none_count = 0
    no_cigarstring_count = 0
    asterisk_for_seq_count = 0
    no_cigarstring_or_asterisk_for_seq_count = 0
    
    supplementary_alignment_count = 0 
    non_autosome_count = 0
    processed_read_count = 0
   
    het_count_list =[]
 
    for read in bam_file: # Each read is an aligned read
        # Count every alignment, before filtering
        total_alignment_count += 1  # only keep if this is being written out somewhere
        
# test read counts
        '''        
        if not read.cigarstring:
            no_cigarstring_count += 1

        if read.seq =="*":
            asterisk_for_seq_count += 1

        if (not read.cigarstring) or read.seq == "*":
            no_cigarstring_or_asterisk_for_seq_count += 1
        ''' 

        # query alignment here means the alignable portion of the read, and therefore excludes clipping, but includes insertions
        # query_alignment_sequence returns None if aligned sequence or cigar string are missing
#        if not read.query_alignment_sequence:
#            query_alignment_none_count += 1
#            continue        
                
        # bitwise flag; 2048 refers to "supplementary alignment"
        if read.flag > 2047:
            supplementary_alignment_count += 1
            continue

        # Skip reads not aligned to "chr{autosome}" 
        if not re.match("chr[0-9]+$", read.reference_name):       
            non_autosome_count += 1
            continue

        processed_read_count += 1
        '''        
    count_file = open("/data/data_repo/castrocp/LRphase/LRphase/Test_output/read_counts.txt", "w")
    count_file.write("total alignments read, before filtering: " + str(total_alignment_count) + "\n")
    count_file.write("total alignments skipped due to having 'None' for 'query_alignment_sequence': " + str(query_alignment_none_count) + "\n")
    count_file.write("not cigarstring: " + str(no_cigarstring_count) + "\n")
    count_file.write("asterisk seq: " + str(asterisk_for_seq_count) + "\n")
    count_file.write("no cigarstring, or asterisk for seq: " + str(no_cigarstring_or_asterisk_for_seq_count) + "\n")
    count_file.write("supplementary alignments filtered: " + str(supplementary_alignment_count) + "\n")
    count_file.write("alignments skipped due to not matching 'chr{autosome}'  pattern: " + str(non_autosome_count) + "\n")
    count_file.write("total alignments moving to processing step (passed all filters): " + str(processed_read_count) + "\n")
        '''    
        # Gap-compressed per-base sequence divergence from minimap2 output
        per_base_error = read.tags["de"][1]       

        # Generate gapped read
        aligned_gapped_read = ''
        gapped_locs = []
        for base, pos in cigar_alignment_padded(read.seq, read.cigarstring, read.pos):
            aligned_gapped_read += base
            gapped_locs.append(pos) # Currently not being used for anything   
        
        # Fetch records from tabix-indexed VCF that overlap with read
        # The "fetch" function takes in start and end as zero-based, half-open coordinates.
        # The very first base of a chromosome is index 0, and the the region includes bases up to, but not including the base at the end coordinate.
        # For example fetch('4', 10, 20) would include all variants overlapping a 10 base pair region from the 11th base of through the 20th base (which is at index 19) of chromosome 4.
        vcf_reader = vcf.Reader(filename=vcf_file)
        het_count = 0
        for vcf_record in vcf_reader.fetch(read.reference_name, read.reference_start, read.reference_end): 
            geno = str(vcf_record.samples[0]["GT"]) # Genotype returned as "1|0", for example.
            # Check for the record being a phased SNP
            if vcf_record.is_snp and vcf_record.heterozygosity == 0.5 and geno[1] == "|":
                snp_pos = vcf_record.POS - read.reference_start - 1  # POS is the one-based position as listed in the VCF file
                read_base = aligned_gapped_read[snp_pos]
                
                # Counter for het sites
                if not read_base == "-":  # don't count sites overlapping deletions in the reference
                    het_count += 1
 
                # Store current SNP info in dictionary
                snp_attributes= {}
                snp_attributes["pos"] = snp_pos
                snp_attributes["gt"] = geno[0] + geno[2]
                snp_attributes["ref"] = str(vcf_record.REF).strip("[").strip("]") # Allele is returned
                snp_attributes["alt"] = str(vcf_record.ALT).strip("[").strip("]")
                snp_attributes["read_base"] = read_base
                snp_attributes["read_ID"] = read.read_name
                snp_attributes["error_rate"] = per_base_error
                
                allele = []
                # Store the DNA base letter corresponding to the VCF notation of ref=0, alt=1
                for i in (0,1):
                    if snp_attributes["gt"][i] == "0":
                        allele.append(snp_attributes["ref"])
                    else:
                        allele.append(snp_attributes["alt"])

                snp_attributes["allele"] = allele
                
                # Add to list of SNP dictionaries
                read_snps.append(snp_attributes)
        
                
        het_count_list.append(het_count)
        
        # Print progress every 10000 reads
        if processed_read_count % 10000 == 0:
            print(f"{processed_read_count} reads have been processed")
    
    het_count_list.sort(reverse = True)
    read_snps_size = len(read_snps)
    unique_readIDs = collections.Counter(d["read_ID"] for d in read_snps)
    num_of_uniquereads = len(collections.Counter(d["read_ID"] for d in read_snps).keys())
    with open(output_dir + "/readID_counts.txt", "a") as counttest:
        for i in het_count_list:    
            counttest.write(str(i) + ' ')
        counttest.write(str(unique_readIDs))
    print(f"Total reads that passed all filters and moved on to VCF step: {processed_read_count}")  #TEST PRINT STATEMENT    
    print(f"Total nubmer of snp_attributes dicts stored in read_snps (total het sites): {read_snps_size}")  #TEST PRINT STATEMENT 
    print(f"Total nubmer of unique read IDs overlapping >= 1 het sites: {num_of_uniquereads}")  #TEST PRINT STATEMENT
    end_process_time=time.time()
    
    total_process_time = end_process_time - start_process_time
    print(f"Finished preparing reads for analysis in {total_process_time:.2f} seconds")

    return(read_snps)
