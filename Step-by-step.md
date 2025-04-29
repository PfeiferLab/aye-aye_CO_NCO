This step-by-step guide provides the context to execute the [`Snakemake` workflows](Snakepit) and seamlessly implement them to other pedigreed variants to identify crossover and non-crossover events.

## Preparing the `Snakemake` environment

[`Snakemake`](https://snakemake.readthedocs.io/en/stable/) is a workflow management system to create reproducible and scalable data analyses. It can be installed with `conda` / `mamba` (along with the required [software](README.md)) as follows:

```bash
mamba create -n recombination_env -c conda-forge -c bioconda snakemake bcftools samtools #different environments / commands can be used if preferred
```

The following files are provided for each of the sections described in the [`README` file](README.md):

- `Snakefile` concatenates the different rules in an output-to-input logic, so the performance of the workflow requires minimal supervision. The parameters and computational resources used for each rule are explicitly stated, so the rules can be easily tweaked and nuanced. The rules invoke software installed in the same `conda` environment or scripts located in the indicated paths. This file should be as interoperable as possible, so all hardcoded information should be included in the `config.yaml`
- `config.yaml` defines the paths, files, software, wildcards and global variables to be used and is customized by the user. It is not required but useful not to hardcode the `Snakefile`.
- `snake_submit.sh` invokes the required files and runs the workflow to the `SLURM` system. The content of this bash file can be seen below.
    - The `SLURM` profile is included in your `Snakemake` installation as indicated [here](https://github.com/Snakemake-Profiles/slurm) ([`cookiecutter`](https://github.com/cookiecutter/cookiecutter) could be useful). Note: newer versions of `Snakemake` (*e.g.* `v.9`) facilitate the `SLURM` integration so these steps may not be necessary.

```bash
snakemake -s Snakefile --configfile config.yaml --profile "slurm" --nolock --rerun-incomplete
```

## Input files

- The [variant calling and filtering worflow](Variant_filtering) requires the following input files:
    - Genotyped `VCF` containing information about `SNPs` and `INDELs` in the pedigree. In [our study](https://academic.oup.com/gbe/advance-article/doi/10.1093/gbe/evaf072/8115344?login=false), the [`Genome Analysis Toolkit`](https://gatk.broadinstitute.org/hc/en-us)'s Best Practice "hard filter" criteria for germline variants has been used. The workflow readily accepts `VCF` files with pedigree information called with other variant callers.
    - Reference genome (`fasta`) of the species of study is required to normalize the `INDELs`
    - Mapped reads (`bam`) are required to calculate the average coverage of each sample
    - A `BED` file including the repetitive regions identified in the `FASTA` file
    - A `PED` pedigree file with the family / trio information
- The [pedigree-based approach](Pedigree_approach) requires the following input file:
    - The filtered `VCF` generated as part of the variant calling and filtering workflow. Importantly, the order of the samples in the `VCF` needs to be mother, father, sample, partner, offspring.
- The [family-based approach](Family_approach) requires the following input file:
    - The filtered `VCF` generated as part of the variant calling and filtering workflow. Importantly, the value of the wildcard offspring (`off`) can be adjusted to any number of available offspring. Our workflow includes `2`, `3` and `4` given our pedigree structure but larger offspring are readily implementable.

## Execution in a nutshell

Provided the input and the `Snakefile` files (`config.yaml` is optional), the only commands the user need to run are:

```bash
mamba activate recombination_env
sh snake_submit.sh
```