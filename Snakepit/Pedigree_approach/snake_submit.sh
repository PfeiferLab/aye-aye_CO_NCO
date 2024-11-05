#!/bin/bash

snakemake -s Snakefile --configfile config.yaml --profile "slurm" --nolock --rerun-incomplete
