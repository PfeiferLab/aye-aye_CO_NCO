## A pedigree-based map of crossovers and non-crossovers in aye-ayes (*Daubentonia madagascariensis*)

Accompanying code for Cyril J. Versoza*, Audald Lloret-Villas*, Jeffrey D. Jensen, Susanne P. Pfeifer. 2025. A pedigree-based map of crossovers and non-crossovers in aye-ayes (*Daubentonia madagascariensis*). [`BioRxiv`](https://www.biorxiv.org/content/10.1101/2024.11.08.622675).

### Code

[`Snakemake`](https://snakemake.readthedocs.io/en/stable/) workflows are provided for the following processes described in the manuscript:

#### Variant filtering

The [variant filtering workflow](Snakepit/Variant_filtering/) defines how the autosomal variants were processed and filtered to obtain a high-confidence `SNP` set. The rules included in this workflow are:

- `raw_indel_VCFs`, retrieves the normalized `INDELs` from the raw genotyped `VCF` with `bcftools norm` and `bcftools view`.
- `bam_coverage`, calculates the sequencing coverage of the `BAM` files with `samtools coverage`.
- `vcf_prep`, concatenates the genotyped autosomal `VCF` files with `bcftools concat` and calls segregating sites with `bcftools view`.
    - The custom script [`vcf_autosomes.sh`](Snakepit/Variant_filtering/scripts/vcf_autosomes.sh) is used to create the list of autosomes required by `bcftools concat`.
- `mask_filter`, provided a `.bed` file with coordinates, `bcftools view` is used to mask the `VCF` file.
- `dp_table`, extracts the read depth (`DP`) in each position and formats the output as a table with `bcftools query`
- `dp_filter`, filters out the variants that have a `DP` less than half or more than twice the average `DP` for that particular sample.
    - The custom script [`DP.sh`](Snakepit/Variant_filtering/scripts/DP.sh) uses `bash` code and `bcftools view` to perform such filtering.
- `gq_filter`, filters out the variants with a genotype quality (`GQ`) lower than `30`, with `bcftools view`.
- `het_filter`, filters out the variants wih an excess of heterozygosity (defined as a `p-value` of `0.01`) with the `bcftools` plug-in `+fill-tags` and `bcftools view`.
- `men_filter`, provded a `.ped` file with pedigree information, it filters out the variants that violate the patterns expected by Mendelian inheritance, with the `bcftools` plugin `+mendelian2`.
- `auto_split`, splits the `VCF` files into autosome-specific `VCF` files, with `bcftools view`.
- `snp_cluster_filter`, removes cluster of variants (defined as `≥ 3 SNPs` within a `10 bp` window)) from the autosomal `VCF` files.
    - The custom script [`SNP_clusters.py`](Snakepit/Variant_filtering/scripts/SNP_clusters.py), developed with `pysam` is used for such purpose.
- `indel_filter`, removes variants located within `10 bp` of an insertion/deletion (`INDEL`).
    - The custom script [`Indels.py`](Snakepit/Variant_filtering/scripts/Indels.py), developed with `pysam` is used for such purpose.
- `chop_ends`, removes variants located within `2 Mb` from the autosome ends.
    - The custom script [`Chop_ends.py`](Snakepit/Variant_filtering/scripts/Chop_ends.py), developed with `pysam` is used for such purpose.

#### Pedigree approach (*under active development*)

- The detection of recombination using the [pedigree approach](Snakepit/Pedigree_approach/)

#### Family approach (*under active development*)

- The detection of recombination events is described using the [family approach](Snakepit/Family_approach/)

#### `Snakemake` execution

Each of the above folder consists of:
- a `Snakefile` with a concatenation of the `Snakemake` rules, the specific parameters and computing resources required
- a `config.yaml` file including the paths and glabl variables, which should be costumized by the user
- a `snake_submit.sh` execution file that triggers the `SLURM` execution of all the rules indicated in the `Snakefile` with the parameters included in the `config.yaml`. Provided the same file names, this execution can be triggered with `snakemake -s Snakefile --configfile config.yaml --profile "slurm" --nolock --rerun-incomplete`.