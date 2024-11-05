# File paths to your input VCF and output VCF
in_vcf = snakemake.input.vcf
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats

import pysam
import subprocess #required for the bash execution

def find_snp_clusters(vcf_file, window_size=10, snp_threshold=3):
    """Finds SNP clusters of more than snp_threshold SNPs within a window_size bp."""
    vcf_in = pysam.VariantFile(vcf_file)  # Open the compressed VCF file (.vcf.gz)
    current_window = []  # List to hold the SNP positions in the current window
    clusters = set()  # Set to hold SNP positions in identified clusters

    # Loop through each record (variant) in the VCF
    for record in vcf_in:
        pos = record.pos
        # Remove SNPs that are outside the 10bp window
        current_window = [p for p in current_window if pos - p <= window_size]
        # Add the current SNP position
        current_window.append(pos)
        # If the number of SNPs in the window exceeds the threshold, mark them as clustered
        if len(current_window) > snp_threshold:
            clusters.update(current_window)

    return clusters

def filter_vcf(vcf_file, output_vcf, clusters):
    """Filters out SNPs in the clusters from the input VCF and writes a new VCF."""
    vcf_in = pysam.VariantFile(vcf_file)  # Open the VCF file
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)  # Open the output VCF

    # Loop through the VCF file and write variants that are not in clusters
    for record in vcf_in:
        if record.pos not in clusters:
            vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

# Step 1: Find clusters of SNPs
snp_clusters = find_snp_clusters(in_vcf)

# Print all positions in the clusters
if snp_clusters:
    print("Positions in identified SNP clusters:")
    for pos in sorted(snp_clusters):
        print(pos)
else:
    print("No SNP clusters found.")

# Step 2: Filter out SNPs in clusters and write to a new VCF file
filter_vcf(in_vcf, out_vcf, snp_clusters)

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