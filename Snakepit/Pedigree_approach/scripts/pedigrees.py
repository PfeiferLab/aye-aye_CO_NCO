input_file = snakemake.input.markers

import pandas as pd

# Define the dataframes and the columns
df = pd.read_csv(input_file, sep='\t', names=["Phase","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Mother","Father","Sample","Partner","Child"])

#Subset the dataframes by dropping some columns and trimming the genotypes
df = df.drop(columns=["ID","QUAL","FILTER","INFO","FORMAT"])
df["Mother"] = df["Mother"].str.split(':',expand=True)[0]
df["Father"] = df["Father"].str.split(':',expand=True)[0]
df["Sample"] = df["Sample"].str.split(':',expand=True)[0]
df["Partner"] = df["Partner"].str.split(':',expand=True)[0]
df["Child"] = df["Child"].str.split(':',expand=True)[0]

df.to_csv(snakemake.output.all_tsv, sep='\t', header=True, index=False) #Writing the results to a tsv file

#df_phase.to_csv(snakemake.output.false_tsv, sep='\t', header=True, index=False) #Writing the results to a tsv file