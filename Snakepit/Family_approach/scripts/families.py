input_file = snakemake.input.markers

import pandas as pd

#different sets of column names based on the number of columns
column_names_2 = ["Scaffold","POS","Reference","Alternate","Maternal","Paternal","Offspring1","Offspring2"]
column_names_3 = ["Scaffold","POS","Reference","Alternate","Maternal","Paternal","Offspring1","Offspring2","Offspring3"]
column_names_4 = ["Scaffold","POS","Reference","Alternate","Maternal","Paternal","Offspring1","Offspring2","Offspring3","Offspring4"]

#function to assign column names...
def import_df(file):
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    #... based on the number of columns
    if len(df.columns) == 8:
        df.columns = column_names_2
    elif len(df.columns) == 9:
        df.columns = column_names_3
    elif len(df.columns) == 10:
        df.columns = column_names_4
    return df

df = import_df(input_file)

#replacing pipes by slashes
df.replace({'0|0': '0/0', '0|1': '0/1', '1|1': '1/1'}, inplace=True) 

#assigning transmitted alleles
condition_mat = ((df['Paternal'] == '0/0') & (df['Offspring1'] == '0/0') | (df['Paternal'] == '1/1') & (df['Offspring1'] == '0/1'))
condition_pat = ((df['Maternal'] == '0/0') & (df['Offspring1'] == '0/0') | (df['Maternal'] == '1/1') & (df['Offspring1'] == '0/1'))

if df.iloc[0]['Maternal'] == '0/1':
    df.loc[condition_mat, 'MatAllOff1'] = 'A'
    df.loc[~condition_mat, 'MatAllOff1'] = 'B'
elif df.iloc[0]['Paternal'] == '0/1':
    df.loc[condition_pat, 'PatAllOff1'] = 'A'
    df.loc[~condition_pat, 'PatAllOff1'] = 'B'

#shared columns to be printed
shared_columns = ['Scaffold', 'POS', 'Reference', 'Alternate', 'Maternal', 'Paternal', 'Offspring1']

#stating whether the other offspring are in phase (SP) or using the other phase (OP) compare to template offspring
for n in range(2, len(df.columns) - 6):
    offspring_column = f"Offspring{n}"
    phase_off_column = f"PhaseOff{n}"
    df.loc[df['Offspring1'] == df[offspring_column], phase_off_column] = 'SP'
    df.loc[df['Offspring1'] != df[offspring_column], phase_off_column] = 'OP'
    InPhaseColumn = f"InPhaseOffspring{n}"
    df.loc[df[phase_off_column] == df[phase_off_column].shift(1), InPhaseColumn] = 'Ph'
    df.loc[df[phase_off_column] != df[phase_off_column].shift(1), InPhaseColumn] = 'PhCh'
    columns_ending_in_n = [col for col in df.columns if col.endswith(str(n))]
    selected_columns = shared_columns + columns_ending_in_n
    df_filtered = df[selected_columns]
    output_file = f"{snakemake.params.offs_tsv}{n}_phased_markers.txt"
    df_filtered.to_csv(output_file, sep='\t', header=True, index=False)

df.to_csv(snakemake.output.all_tsv, sep='\t', header=True, index=False)