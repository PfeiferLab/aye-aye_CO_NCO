markers_file = snakemake.input.markers

import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

# Reading the phase-informative markers and import them as a dataframe
df = pd.read_csv(markers_file, sep='\t')

# Drop last column and re-name the Phase column
df = df.drop(df.filter(like='InPhaseOffspring').columns, axis=1)
df = df.rename(columns={col: 'Phase' for col in df.columns if col.startswith('PhaseOff')})

# Create a column to detect changes in haplotype phase
df['InPhaseAfter'] = df.Phase.eq(df.Phase.shift(-1))
if not df.empty:
    df.loc[df.index[-1], 'InPhaseAfter'] = True

# Define the distance threshold and the minimum number of points
threshold = 5000
min_points = 4

# Function to remove rows within a window
def remove_close_rows(df, threshold=5000, min_points=4):
    filtered_df = df[df['InPhaseAfter'] == False]

    # Initialize a list to store indices to be removed from the original df
    indices_to_remove = []

    # Iterate through the filtered dataframe
    for i in range(len(filtered_df)):
        current_position = filtered_df.iloc[i]['POS']
        
        window = filtered_df[(filtered_df['POS'] >= current_position) & (filtered_df['POS'] <= current_position + threshold)]

        if len(window) >= min_points:
            indices_to_remove.extend(df[(df['POS'] >= current_position) & (df['POS'] <= current_position + threshold)].index)
    
    # Capture the rows to be removed and remove them
    deleted_rows = df.loc[indices_to_remove]
    df_cleaned = df.drop(index=indices_to_remove)

    return df_cleaned, deleted_rows

# Apply the function to clean the dataframe and get deleted rows
df_clean, deleted_rows = remove_close_rows(df)

# Reset the index after filtering and write the results to tsv files
df_clean.reset_index(drop=True, inplace=True)
#deleted_rows.to_csv(snakemake.output.deleted_markers, sep='\t', header=True, index=False)
df_clean.to_csv(snakemake.output.clean_markers, sep='\t', header=True, index=False)

# Extra columns are removed
df_phase = df_clean.drop(columns=['InPhaseAfter']).copy()

# Creating two new columns to check changes of phase
df_phase['InPhaseBefore'] = df_phase.Phase.eq(df_phase.Phase.shift(1))
if not df_phase.empty: 
    df_phase.loc[df_phase.index[0], 'InPhaseBefore'] = True
df_phase['InPhaseAfter'] = df_phase.Phase.eq(df_phase.Phase.shift(-1))
if not df_phase.empty: 
    df_phase.loc[df_phase.index[-1], 'InPhaseAfter'] = True 

# Breakpoints (changes of phase) are selected
df_bkpoint = df_phase[(df_phase["InPhaseBefore"] == False) | (df_phase["InPhaseAfter"] == False)]
df_SNP = df_phase[(df_phase['InPhaseBefore'] == False) & (df_phase['InPhaseAfter'] == False)]
df_SNP.to_csv(snakemake.output.snps, sep='\t', header=True, index=False)

# A dataframe of the breakpoints encompassing more than 1 SNP is created
df_ChPh = df_bkpoint[~df_bkpoint.POS.isin(df_SNP.POS)]
#df_ChPh.to_csv(snakemake.output.ChPh_del, sep='\t', header=True, index=False)

# New columns
df_ChPh_trim = df_ChPh.drop(columns=['InPhaseBefore','InPhaseAfter']).copy()

df_ChPh_trim['InPhaseBefore'] = df_ChPh_trim.Phase.eq(df_ChPh_trim.Phase.shift(1))
if not df_ChPh_trim.empty: 
    df_ChPh_trim.loc[df_ChPh_trim.index[0], 'InPhaseBefore'] = True 
df_ChPh_trim['InPhaseAfter'] = df_ChPh_trim.Phase.eq(df_ChPh_trim.Phase.shift(-1))
if not df_ChPh_trim.empty: 
    df_ChPh_trim.loc[df_ChPh_trim.index[-1], 'InPhaseAfter'] = True

df_ChPh_trim = df_ChPh_trim[(df_ChPh_trim["InPhaseBefore"] == False) | (df_ChPh_trim["InPhaseAfter"] == False)]

df_MidPoint = df_ChPh_trim.copy()
df_MidPoint = df_MidPoint.reset_index(drop=True)

# A "ResMidPoint" column is created and the values are propagated to both bounds of the Midpoint
df_MidPoint['ResMidPoint'] = df_MidPoint.apply(
    lambda row: (row['POS'] + df_MidPoint.loc[row.name - 1, 'POS']) / 2
    if not row['InPhaseBefore'] and row['InPhaseAfter'] and row.name > 0 else None, axis=1
)
df_MidPoint['ResMidPoint'] = df_MidPoint['ResMidPoint'].bfill()

df_MidPoint.to_csv(snakemake.output.ChPh, sep='\t', header=True, index=False) # Writing the results to a tsv file