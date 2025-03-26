in_vcf = snakemake.input.vcf
out_vcf = snakemake.output.vcf
stats_vcf = snakemake.output.stats

import pysam
import subprocess #required for the bash execution

def find_snp_clusters(vcf_file, window_size=10, snp_threshold=3):
    """Finds SNP clusters of more than snp_threshold SNPs within a window_size bp."""
    vcf_in = pysam.VariantFile(vcf_file)
    current_window = []
    clusters = set()

    for record in vcf_in:
        pos = record.pos
        current_window = [p for p in current_window if pos - p <= window_size]
        current_window.append(pos)
        #when the number of SNPs in the window exceeds the threshold, mark as clustered
        if len(current_window) > snp_threshold:
            clusters.update(current_window)

    return clusters

def filter_vcf(vcf_file, output_vcf, clusters):
    """Filters out SNPs in the clusters from the input VCF and writes a new VCF."""
    vcf_in = pysam.VariantFile(vcf_file)
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    #write variants that are not in clusters
    for record in vcf_in:
        if record.pos not in clusters:
            vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

#find SNP clusters
snp_clusters = find_snp_clusters(in_vcf)

#print positions in clusters
if snp_clusters:
    print("Positions in identified SNP clusters:")
    for pos in sorted(snp_clusters):
        print(pos)
else:
    print("No SNP clusters found.")

#filter out SNPs in clusters
filter_vcf(in_vcf, out_vcf, snp_clusters)

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