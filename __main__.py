"""Main Script for LRphase"""

import bamnostic as bm
import os
import sh
from sh import ErrorReturnCode
import click
from data_compile import process_data
from data_analysis import count_alleles
import pysam
import time

def run_alignment(ref, reads, output_dir):
    """Align specified reads file to reference genome via minimap2."""
    try:
        cmd = sh.Command("./exec/minimap2/minimap2")
        cmd("-ax", "map-ont", "-L", "--secondary=no", "--sam-hit-only", ref, reads, _out=output_dir + "/alignment.sam") 
        # "-a" = creates aligned SAM file, "x" to choose a preset (map-ont in this case) 
        # "-L" option is used when working with ultra-long nanopore reads to account for CIGAR strings > 65,535 characters
        # map-ont uses ordinary minimizers as seeds and is recommended for use with Oxford Nanopore reads
    
    except ErrorReturnCode:
        print("Error: ", ErrorReturnCode)
        raise Exception("Failed to align genome via Minimap2... Exiting")
    except KeyboardInterrupt:
        kill = sh.Command("pkill")
        kill.bake("-lf")
        kill("minimap2")
        raise Exception("Manual Abort From User, \
                         Minimap2 has stopped... Exiting")

    # Convert SAM to BAM and index BAM file
    pysam.view("-bS" , output_dir + "/alignment.sam", "-o", output_dir + "/alignment.bam", catch_stdout=False)
    bam_file = output_dir + "/alignment.bam"
    pysam.sort(bam_file, "-o",  output_dir + "/alignment.sorted.bam", catch_stdout=False)
    sorted_bam = output_dir + "/alignment.sorted.bam"
    pysam.index(sorted_bam)

    
# Configure command line options
@click.command()
@click.option(
    "--output", "-o", "output_directory", #default="output_insertion",
    help="Output directory", #default=output",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
)
@click.option(
    "--ref", "-r", "reference_genome", #default="./data/reference_hg38.fna",
    help="Reference Genome File", #default=./data/reference_hg38.fna",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)
@click.option(
    "--fastq", "-f", "long_read_fastq", #default="./data/1753.fastq",
    help="Long Read FastQ File", #default=./data/1753.fastq",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)
@click.option(
    "--vcf", "-v", "vcf_file", #default="./data/hg001.vcf",
    help="VCF File with Tabix index", #default=./data/hg001.vcf",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)


def main(output_directory, reference_genome, long_read_fastq, vcf_file):

    """Create output subdirectories."""
    print("************* Creating output subdirectories")
    source_dir = os.path.dirname(os.path.realpath(__file__))
    output_dir = output_directory 

    alignments_dir = os.path.join(output_directory, "alignment-output")
    if os.path.exists(alignments_dir):
        print("Alignments Subdirectory Already Exists")
    else:
        os.makedirs(alignments_dir)
    
    results_dir = os.path.join(output_directory, "results")
    if os.path.exists(results_dir):
        print("Results Subdirectory Already Exists")
    else:
        os.makedirs(results_dir)

    alignment_start_time = time.time()
    # Can be set to "False" if alignement files already exist in output_directory/alignment-output
    align = True
    print("************* Beginning Sequence Alignment")
    if align:
        run_alignment(reference_genome, long_read_fastq, alignments_dir)
    
    alignment_end_time = time.time()
    time_alignment = alignment_end_time-alignment_start_time
    print(f"Alignment and conversion of SAM to BAM finished in {time_alignment:.2f} seconds")

    read_snps = process_data(output_dir, reference_genome, vcf_file)

    count_alleles(read_snps, results_dir)


if __name__ == "__main__":
    main()
