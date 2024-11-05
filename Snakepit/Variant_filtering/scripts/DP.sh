#!/usr/bin/env bash

pat='autosomes_mask_DP_sample'

cp ${snakemake_input[vcf]} ${snakemake_output[vcf_pre]}; cp ${snakemake_input[tbi]} ${snakemake_output[tbi_pre]}; 

samples=${snakemake_params[samples]}
for ((sam=1; sam<=$samples; sam++))
do
	col=$(expr ${sam} + 2)
	pycol=$(expr ${sam} - 1)
	cov=`awk -vcol="$col" '{c+=$col;++n} END {print c/n}' < ${snakemake_input[table]}`
	maxcov=$(bc <<< "$cov * 2")
	mincov=$(bc <<< "$cov * 0.5")
	bcftools view ${snakemake_output[path]}/${pat}_${pycol}.vcf.gz -i "FORMAT/DP[$pycol] > $mincov && FORMAT/DP[$pycol] < $maxcov" -Oz -o ${snakemake_output[path]}/${pat}_${sam}.vcf.gz
	tabix -fp vcf ${snakemake_output[path]}/${pat}_${sam}.vcf.gz
	bcftools stats ${snakemake_output[path]}/${pat}_${sam}.vcf.gz > ${snakemake_output[path]}/${pat}_${sam}.stats
	echo 'Sample:' $sam' // Coverage:' $cov' // Maximum coverage:' $maxcov' // Minimum coverage:' $mincov' // Column:' $col' // PyColumn:' $pycol
	echo '--------------'
done

cp ${snakemake_output[vcf_post]} ${snakemake_output[vcf]}; cp ${snakemake_output[tbi_post]} ${snakemake_output[tbi]}; cp ${snakemake_output[stats_post]} ${snakemake_output[stats]}
