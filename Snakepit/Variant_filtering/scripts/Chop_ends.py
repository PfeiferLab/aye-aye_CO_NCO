in_vcf = snakemake.input.vcf
index_file = snakemake.input.fai
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats
chrom = snakemake.wildcards.chrom

import pysam
import subprocess # Required for the bash execution

def get_upper_limit(index_file_path, chrom): #autosome length are obtained from the fasta file
    upper_limit = None
    with open(index_file_path, 'r') as f:
        for line in f:
            if line.startswith(f"prefix_{chrom}"):
                columns = line.strip().split('\t')
                if len(columns) > 1:
                    upper_limit = int(columns[1]) - 2_000_000
                break 
    return upper_limit

def filter_vcf(input_vcf_path, output_vcf_path, upper_limit):
    with pysam.VariantFile(input_vcf_path) as vcf_in:
        #output VCF file with the same header
        with pysam.VariantFile(output_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                position = record.pos
                
                #check if position is within range
                if 2_000_000 < position < upper_limit:
                    vcf_out.write(record)

upper_limit = get_upper_limit(index_file, chrom)
filter_vcf(in_vcf, out_vcf, upper_limit)

def run_tabix(out_vcf):
    """Index a VCF file."""
    command = ['tabix', '-fp', 'vcf', out_vcf]

    result = subprocess.run(command, capture_output=True, text=True)

def run_bcftools_stats(out_vcf, stats_vcf):
    """Stats of the VCF file"""
    command = ['bcftools', 'stats', out_vcf]

    with open(stats_vcf, 'w') as f:
        result = subprocess.run(command, stdout=f, stderr=subprocess.PIPE, text=True)

#run tabix and bcftools stats
run_tabix(out_vcf)
run_bcftools_stats(out_vcf, stats_vcf)