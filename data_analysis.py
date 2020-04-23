"""Data analysis functions for LRphase"""

from operator import itemgetter
import pandas as pd
import numpy as np
import math
import time

def calc_likelihood_ratio(error, pat_matches, mat_matches):
    p_hpat = (1-error)**pat_matches * (error/3)**mat_matches
    if p_hpat == 0:
        p_hpat = 1.0e-322 # stop calculations as the value approaches zero
    
    p_hmat = (1-error)**mat_matches * (error/3)**pat_matches
    if p_hmat == 0:
        p_hmat = 1.0e-322
    
    l_ratio = p_hpat/p_hmat
    if l_ratio == 1/1e-309:
        l_ratio = 1e+308

    return( math.log10(l_ratio) )

def create_output(reads, hets, results_dir, het_locs, error_rates):
    """Create output file"""
    print("************* Phasing and writing results file")
    
    # This will store all the counts/statistics for each read, to be written to file
    all_rows = []

    for read in reads:
        minimum_error_rate = error_rates[read]
        mat = 0
        pat = 0
        het = hets[read] #total number of het sites overlapping this read

        for snp, loc in zip(reads[read], het_locs[read]):
            # Count the number of alleles matching the paternal/maternal reference haplotype
            if snp == 1:
                pat = pat + 1
            if snp == 2:
                mat = mat + 1
        
        # Calculate updated likelihood ratio
        likelihood_ratio = calc_likelihood_ratio(minimum_error_rate, pat, mat)

        # Store counts/statistics for this read
        row = []
        row.append(read) # read name
        row.append(mat) # number of SNPs matching maternal reference haplotype
        row.append(pat) 
        row.append(het) # total number of heterozygous sites within the read
        row.append(minimum_error_rate)
        row.append(likelihood_ratio)
        all_rows.append(row)
 
    # sorts reads by "key" which is the 4th item, called by "itemgetter(3)".  In this case, total het sites.
    all_rows = sorted(all_rows, key=itemgetter(3), reverse=True)
    if len(all_rows) == 0:
        all_rows = [["none", "none", "none", "none", "none", "none"]]
    
    # Write final output file
    df = pd.DataFrame(np.array(all_rows),
                       columns=['read name', 'maternal matches',
                                'paternal matches', 'heterozygous positions',
                                'Empirical Error Rate', 'log10_likelihood_ratio'])

    with open(results_dir + "/phasing_stats.tsv", "w") as tsv_file:
        df.to_csv(path_or_buf=tsv_file, sep="\t", index=False)


def count_alleles(read_snps, results_dir):
    start_allele_count_time = time.time()
    print("************* Begin counting alleles")
   
    reads = {}
    hets = {}
    het_locs = {}
    error_rates = {}
 
    # read_snps is a list of dictonaries
    # Each dictionary contains the information for one SNP
    for snp_info in read_snps:
        if snp_info["read_ID"] not in reads:
            reads[snp_info["read_ID"]] = []
            hets[snp_info["read_ID"]] = 0
            het_locs[snp_info["read_ID"]] = []
            error_rates[snp_info["read_ID"]] = snp_info["error_rate"]        
    
        # Count the number of heterozygous sites
        # genotype at each read is "0" for reference allele, "1" for alternate allele
        if snp_info["gt"] == "01" or snp_info["gt"] == "10":
            if not snp_info["read_base"] == "-": # ignore positions overlapping deletions
                hets[snp_info["read_ID"]] = hets[snp_info["read_ID"]] + 1
                het_locs[snp_info["read_ID"]].append(snp_info["pos"])
        
        # Check that the read base is one of the alleles listed in the VCF record at that site
        if snp_info["read_base"] in snp_info["allele"]:
            
            # Check for the different genotype scenarios 
            if snp_info["gt"] == "01": #ref/alt
                if snp_info["read_base"] == snp_info["alt"]:
                    # If base == alternate allele, the SNP is the 2nd allele
                    reads[snp_info["read_ID"]].append(2)
                elif snp_info["read_base"] == snp_info["ref"]:
                    reads[snp_info["read_ID"]].append(1)

            elif snp_info["gt"] == "10":  # alt/ref
                if snp_info["read_base"] == snp_info["alt"]:
                    # if base == alt, SNP is the 1st allele
                    reads[snp_info["read_ID"]].append(1)

                elif snp_info["read_base"] == snp_info["ref"]:
                    reads[snp_info["read_ID"]].append(2)
            
        # count the het site even if the read allele is different from the ones in the VCF
        else:
            if not snp_info["read_base"] == "-":
                hets[snp_info["read_ID"]] = hets[snp_info["read_ID"]] + 1


    end_allele_count_time = time.time()
    total_allele_count_time = end_allele_count_time - start_allele_count_time
    print(f"Finished counting alleles in {total_allele_count_time:.2f} seconds")

    create_output(reads, hets, results_dir, het_locs, error_rates)
