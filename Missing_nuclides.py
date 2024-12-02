"""
Missing_nuclides.py

This code checks to see if any of the nuclides in the ENDF database are not in the ORION database.
This is done by reading in the .csv file that contains the ENDF nuclides and the .xlsx file that contains the ORION nuclides. These are loaded as pandas dataframes.
The pandas dataframes are compared and saved as new dataframe, which contains all of the ENDF nuclides in the 1st column and the ORION nuclides in the 2nd column. If the nuclide is not in the ORION database then the element in the 2nd column is 'Not in ORION'.
The final dataframe is outputted to a csv file in the working directory.

Author: Sam Taylor (sam.taylor@newcleo.com)
Last Edited: 22/11/2024
"""
# ---------  Import Libraries  --------- #
import numpy as np
import pandas as pd
# -------------  Main Code  ------------ #

# Read in the .csv file containing the ENDF nuclide information and stores it as a pandas dataframe. Assigns the first column to the nuclide name variable
dir_to_endf_nuclides = r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/Generating_MPR_file/LFR30_Reaction_Data/ZAID_results.csv'
df_endf = pd.read_csv(dir_to_endf_nuclides)
nuclides_in_endf = df_endf[df_endf.columns[0]]

# Read in the .xlsx file containing the ORION nuclide information and stores it as a pandas dataframe.
dir_to_orion_nuclides =r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/Generating_MPR_file/orion_nuclides_list.xlsx'
df_ORION_nuclides = pd.read_excel(dir_to_orion_nuclides, header = None)
# Columns in excel file do not contain headers so create them here
df_ORION_nuclides.columns = ['Nuclide Name', 'Buffer Mass']

# Remove the N- prefix from the names of the nuclides in the ORION database. This converts column 0 to an array and so afterwards convert back to pandas dataframe
df_nuclide_without_prefix = pd.DataFrame({'ORION Nuclides': [element.split('-', 1)[1] for element in df_ORION_nuclides['Nuclide Name']]})

# Create a new dataframe with 2 columns. 1st column is the ENDF nuclides and the 2nd column is for the ORION nuclides, which default to 'Not in ORION'
comparison_df = pd.DataFrame({'ENDF Nuclides': nuclides_in_endf})
comparison_df['ORION Nuclides'] = 'Not in ORION'

# Checks to see if the ORION nuclides are in the 1st column of comparison_df and appends them to a new dataframe if they are
matching_nuclides = df_nuclide_without_prefix[df_nuclide_without_prefix['ORION Nuclides'].isin(comparison_df['ENDF Nuclides'])]
# Iterates through each row in comparison_df and appends the matching ORION nuclides to the correct row.
for index, row in comparison_df.iterrows():
    if row['ENDF Nuclides'] in matching_nuclides['ORION Nuclides'].values:
        comparison_df.loc[index, 'ORION Nuclides'] = row['ENDF Nuclides']

# Save the whole dataframe to a csv file
file_name = 'Missing_Nuclides.csv'
comparison_df.to_csv(file_name, index = False)