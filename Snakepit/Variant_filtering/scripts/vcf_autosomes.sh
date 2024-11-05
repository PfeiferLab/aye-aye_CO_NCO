vcf_path='/path/to/file_folder'
SNP_path='/path/to/file_folder'

for auto in {1..8} {10..15} 
do
	echo ${vcf_path}/Chr_${auto}.vcf >> ${SNP_path}/autos.txt
done
