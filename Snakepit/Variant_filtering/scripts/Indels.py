snp_vcf = snakemake.input.snp_vcf
indel_vcf = snakemake.input.indel_vcf
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats

import pysam
import subprocess  # Required for the bash execution

def nearby_indels(snp_vcf_path, indel_vcf_path, output_vcf_path, distance=10):
    snp_vcf = pysam.VariantFile(snp_vcf_path)
    indel_vcf = pysam.VariantFile(indel_vcf_path)

    with pysam.VariantFile(output_vcf_path, 'w', header=snp_vcf.header) as output_vcf:
        #list of INDELs
        indel_positions = []
        for indel_record in indel_vcf:
            indel_positions.append(indel_record.pos)

        indel_set = set(indel_positions)

        for snp_record in snp_vcf:
            snp_pos = snp_record.pos
            
            #check INDELs within the specified distance
            snp_in_range = False
            for offset in range(-distance, distance + 1):
                if (snp_pos + offset) in indel_set:
                    snp_in_range = True
                    print(f"SNP at position {snp_pos} is within range of an INDEL at position {snp_pos + offset}")
                    break
            
            # SNPs to a new VCF if no nearby INDEL
            if not snp_in_range:
                output_vcf.write(snp_record)

nearby_indels(snp_vcf, indel_vcf, out_vcf, distance=10)

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