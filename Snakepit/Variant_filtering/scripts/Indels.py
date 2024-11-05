# File paths to your input VCF and output VCF
snp_vcf = snakemake.input.snp_vcf
indel_vcf = snakemake.input.indel_vcf
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats

import pysam
import subprocess  # Required for the bash execution

def nearby_indels(snp_vcf_path, indel_vcf_path, output_vcf_path, distance=10):
    # Open SNP and INDEL VCF files with pysam
    snp_vcf = pysam.VariantFile(snp_vcf_path)
    indel_vcf = pysam.VariantFile(indel_vcf_path)

    # Create an output VCF file
    with pysam.VariantFile(output_vcf_path, 'w', header=snp_vcf.header) as output_vcf:
        # Create a list of INDEL positions for quick access
        indel_positions = []
        for indel_record in indel_vcf:
            indel_positions.append(indel_record.pos)

        # Convert to set for faster lookups
        indel_set = set(indel_positions)

        # Iterate through SNPs in the SNP VCF file
        for snp_record in snp_vcf:
            snp_pos = snp_record.pos
            
            # Check for INDELs within the specified distance
            snp_in_range = False
            for offset in range(-distance, distance + 1):
                if (snp_pos + offset) in indel_set:
                    snp_in_range = True
                    print(f"SNP at position {snp_pos} is within range of an INDEL at position {snp_pos + offset}")
                    break
            
            # Write the SNP to the new VCF if no nearby INDEL found
            if not snp_in_range:
                output_vcf.write(snp_record)

# Run the code
nearby_indels(snp_vcf, indel_vcf, out_vcf, distance=10)

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