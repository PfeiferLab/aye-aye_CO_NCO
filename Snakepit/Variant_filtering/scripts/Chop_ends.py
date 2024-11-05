# File paths to your input VCF and output VCF
in_vcf = snakemake.input.vcf
index_file = snakemake.input.fai
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats
chrom = snakemake.wildcards.chrom

import pysam
import subprocess # Required for the bash execution

def get_upper_limit(index_file_path, chrom):
    upper_limit = None
    with open(index_file_path, 'r') as f:
        for line in f:
            if line.startswith(f"prefix_{chrom}"):
                columns = line.strip().split('\t')  # Split by tabs
                if len(columns) > 1:
                    upper_limit = int(columns[1]) - 2_000_000 # Get the second column and subtract 2,000,000
                break  # Stop searching once we find the relevant line
    return upper_limit

def filter_vcf(input_vcf_path, output_vcf_path, upper_limit):
    # Open the input VCF file with pysam
    with pysam.VariantFile(input_vcf_path) as vcf_in:
        # Create an output VCF file with the same header
        with pysam.VariantFile(output_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            # Iterate through records in the input VCF
            for record in vcf_in:
                position = record.pos
                
                # Check if the position is within the specified range
                if 2_000_000 < position < upper_limit:
                    # Write the record to the output VCF file if within range
                    vcf_out.write(record)

# Get the upper limit from the third file
upper_limit = get_upper_limit(index_file, chrom)
filter_vcf(in_vcf, out_vcf, upper_limit)

def run_tabix(out_vcf):
    """Runs the tabix command to index a VCF file."""
    command = ['tabix', '-fp', 'vcf', out_vcf]  # Create the command as a list

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error running tabix: {result.stderr}")
    else:
        print(f"Successfully indexed {out_vcf} with tabix.")

def run_bcftools_stats(out_vcf, stats_vcf):
    """Runs the bcftools stats command on a VCF file."""
    command = ['bcftools', 'stats', out_vcf]  # Create the command as a list

    # Open the stats_vcf file in write mode
    with open(stats_vcf, 'w') as f:
        # Run the command and redirect output to stats_vcf
        result = subprocess.run(command, stdout=f, stderr=subprocess.PIPE, text=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools stats: {result.stderr}")
    else:
        print(f"Successfully generated stats for {out_vcf} in {stats_vcf}.")

# Run tabix to index the VCF file
run_tabix(out_vcf)

# Run bcftools stats to generate statistics for the VCF file
run_bcftools_stats(out_vcf, stats_vcf)