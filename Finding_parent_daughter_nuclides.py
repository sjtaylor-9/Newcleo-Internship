"""
Finding_parent_daughter_nuclides.py

This code is responsible for finding the ZAIDs of the parent and daughter nuclides, which are needed for the MPR file. The ZAIDs of the parent and daughter nuclides for each reaction are outputted to a .csv file.
In addition, the ORION IDs NumberParents, and NumberReactions are also written to the .csv file. These are also all neded for the MPR file.
The code calculates the parent/daughter ZAIDs for each nuclide and then determines if these exist in both the ENDF and ORION databases. All of the reaction information for the nuclide is then appended to a dictionary and written to a .csv file using pandas, once all nuclides that exist in both the ENDF and ORION databases have been examined.

Author: Sam Taylor (sam.taylor@newcleo.com)
Last Edited: 22/11/2024
"""
# ---------  Import Libraries  --------- #
import numpy as np
import pandas as pd
import argparse
import os
# -------------  Functions  ------------ #
def parse_arguments():
    """
    Parses the arguments needed along the code. Arguments:

    --path                  Used to specify the directory in which the one-group cross-section data should be written. It is not required,
                            in the case it is not specified, the default path is the current working directory.
    --cross_section_data    Used to specify the directory in which the public data regarding the ENDF and JEFF reaction cross-sections can be found.
    --reaction_data     
    
    Returns the parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reactor",
        type = str,
        choices = ['LFR30', 'LFR200'],
        required = True,
        help = 'flag to set whether the LFR30 or LFR200 is to be examined'
    )
    parser.add_argument(
        "--enrichment_zone",
        type = str,
        choices = ['inner', 'middle', 'outer', 'FA'],
        required = True,
        help = 'flag to set which enrichment zone in the LFR200 is to be examined'
    )
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--cross_section_data",
        type=dir_path,
        required=True,
        help="flag to set the path where the input cross-section data should be found"
    )
    parser.add_argument(
        "--reaction_data",
        type=dir_path,
        required=True,
        help="flag to set the path where the input reaction data should be found"
    )
    parser.add_argument(
        "--burnup",
        type = str,
        choices = ['BoL', 'EoL', 'single'],
        required = True,
        help = 'flag to set which enrichment zone in the LFR200 is to be examined'
    )
    return parser.parse_args()

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def build_nuclide_name(nuclide):
    """
    Constructs the name of the nuclide in the form PU243 from it's ZAID.
    The ZAID is in the form ZZAAA, so the function sets Z to the first 2 digits and A as the rest, whilst removing all leading/trailing zeroes.

    Args:
        nuclide_string (str): The ZAID number of the nuclide in the form ZZAAA.

    Returns:
        nuclear_notation (str): The name of the parent/daughter nuclide in the form PU243.
        atomic_number (int): The atomic number of the parent/daughter nuclide.
        mass_number (int): The mass number of the parent/daughter nuclide.
    """
    nuclide = int(nuclide)
    # Atomic number can be extracted from the ZAID by dividing by 10,000 and ignoring the remainder
    atomic_number = nuclide // 10000
    # Mass number can be extracted from the ZAID by taking the remainder when dividng 10,000 and then dividing this by 10
    mass_number = nuclide % 10000
    mass_number = mass_number // 10

    # Nuclide name in form PU243
    nuclear_notation = element_symbols[int(atomic_number)-1].upper() + str(mass_number)

    if (nuclide % 10 == 1):
        nuclear_notation = nuclear_notation + "(1)"
    elif (nuclide % 10 == 2):
        nuclear_notation = nuclear_notation + "(2)"
    
    return nuclear_notation, atomic_number, mass_number

def build_nuclide_ID(atomic_number, mass_number):
    """
    Constructs the ZAID of the nuclide from it's atomic and mass numbers. The ZAIDs are used to identify the nuclides.
    The ZAID is in the form ZZAAAI, where Z is the atomic number and A is the mass number. I > 0 identifies that an isotope is metastable.

    Args:
        atomic_number (int): The atomic number of the nuclide, Z.
        mass_number (int): The mass number of the nuclide, A.

    Returns:
        ID (int): The ZAID number of the nuclide in the form ZZAAAI (I = 0 is the ground state nuclide and I > 0 for metastable isomers).
    """
    # Nuclear IDs in form ZZAAAI, so calculated as 10,000Z + A * 10 + I, where I = 0 is the ground state and I = 1 and 2 are metastable isomers.
    ID = 1000 * int(atomic_number)
    ID = ID + int(mass_number)
    ID = ID * 10
    
    return ID

def calculate_parent(atomic_number, mass_number, nuclide_in_dataframe):
    """
    Determines if the parent nuclide of a given reaction exists in both the ENDF and ORION databases and if so then the ZAID of the parent nuclide is calculated.
    If the parent nuclide does not exist in either databse then the ZAID is set to 'NA'.
    The parent of a fission product (reaction 18) is too difficult to determine and so is not considered.
    
    Args:
        atomic_number (int): The atomic number (Z) of the nuclide under examination (i.e. the daughter of the possible parent nuclide).
        mass_number (int): The mass number (A) of the nuclide under examination (i.e. the daughter of the possible parent nuclide).
        nuclide_in_dataframe (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.

    Returns:
        parent_ID_MT (int or str): Where (MT = 16, 17, 102, 103). If the parent nuclide of the specified reaction exists in both the ENDF and ORION databases, then the ZAID is outputted. However, if it does not then the ZAID is set to 'NA'. 
        NumberParents (int): The total number of parent nuclides that can undergo reactions any of 16, 17, 102 and 103 to create the nuclide.
    """
    # If nuclide is daughter of reaction 16 (n,2n) then parent should have A+1
    parent_A_16 = mass_number + 1
    parent_16 = build_nuclide_ID(atomic_number, parent_A_16)
    
    # If nuclide is daughter of reaction 17 (n,3n) then parent should have A+2
    parent_A_17 = mass_number + 2
    parent_17 = build_nuclide_ID(atomic_number, parent_A_17)

    # If nuclide is daughter of reaction 102 then parent should have A-1
    parent_A_102 = mass_number - 1
    parent_102 = build_nuclide_ID(atomic_number, parent_A_102)

    # If nuclide is daughter of reaction 103 then parent should have same A and Z+1
    parent_Z_103 = atomic_number + 1
    parent_103 = build_nuclide_ID(parent_Z_103, mass_number)

    # Checks to see if the parent nuclide is in the database by comparing the possible parent ID with the ones in the first column of the .csv file
    NumberParents = 0
    NumberParents, parent_16_exists = find_parent(NumberParents, parent_16, nuclide_in_dataframe, MT_value = 16)
    if parent_16_exists == False:
        parent_16 = 'NA'

    NumberParents, parent_17_exists = find_parent(NumberParents, parent_17, nuclide_in_dataframe, MT_value = 17)
    if parent_17_exists == False:
        parent_17 = 'NA'
    
    NumberParents, parent_102_exists = find_parent(NumberParents, parent_102, nuclide_in_dataframe, MT_value = 102)
    if parent_102_exists == False:
        parent_102 = 'NA'
    
    NumberParents, parent_103_exists = find_parent(NumberParents, parent_103, nuclide_in_dataframe, MT_value = 103)
    if parent_103_exists == False:
        parent_103 = 'NA'

    return NumberParents, parent_16, parent_17, parent_102, parent_103

def find_parent(NumberParents, parent_id, df_containing_nuclides, MT_value):
    """
    Determines if the parent nuclide exists in both the ENDF and ORION databases.
    This is done by taking in the possible parent ZAID as calculated from the calculate_parent() function and comparing it to the list of all the reactions in the ENDF database.

    Args:
        NumberParents (int): The running total of the number of parents that create the nuclide.
        parent_id (int): The ZAID of the possible parent nuclide.
        df_containing_nuclides (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.
        MT_value (int): The MT value of the given reaction (16, 17, 18, 102, 103).

    Returns:
        NumberParents (int): The running total of the number of parents that create the nuclide.
        parent_exists (bool): Boolean operator that is defaulted to False. Set to True only if the parent nuclide exists in both the ENDF and ORION databases.
    """
    parent_exists = False
    # Iterates through all the rows of the array containing the ENDF reaction data
    for row in reaction_data:
        # If the row of the array has the same ZAID as the parent nuclide and the same reaction number then the parent exists in the ENDF database
        if row[0] == parent_id and row[1] == MT_value:
            parent_exists = True
            # Outputs the boolean operator to determine if the parent nuclide exists in the ORION database (defaults to False)
            parent_exists = if_parent_daughter_exist_in_orion(parent_id, df_containing_nuclides, parent_exists)
            if parent_exists == True:
                # The number of parents total increases by 1 if the parent exists in the ORION database and the loop is exited
                NumberParents = NumberParents + 1
                break
            
    return NumberParents, parent_exists

def find_daughter(daughter_id, df_containing_nuclides):
    """
    Determines if the daughter nuclide exists in both the ENDF and ORION databases.
    This is done by taking in the possible parent ZAID as calculated from the calculate_daughter() function and comparing it to the list of all the reactions in the ENDF database.
    The MT value of the reaction is not needed since the function is only called when the reaction is found in the ENDF database.

    Args:
        daughter_id (int): The ZAID of the possible daughter nuclide.
        df_containing_nuclides (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.

    Returns:
        daughter_exists (bool): Boolean operator that is defaulted to False. Set to True only if the daughter nuclide exists in both the ENDF and ORION databases
    """
    daughter_exists = False

    # Iterates through all the rows of the array containing the ENDF reaction data
    for row in reaction_data:
        # If the row of the array has the same ZAID as the daughter nuclide then the daughter exists in the ENDF database
        if row[2] == daughter_id:
            daughter_exists = True
            # Checks to see if the daughter nuclide exists in the ORION database
            daughter_exists = if_parent_daughter_exist_in_orion(daughter_id, df_containing_nuclides, daughter_exists)
            break
    
    return daughter_exists

def find_ORION_ID(name_of_nuclide, df_containing_nuclides, nuclide_in_ORION):
    """
    Determines if the examined nuclide exists in the ORION database.
    This is done because the U241 nuclide exists in the ENDF database but not the ORION one. This is the only case of this but the process has been generalised for validation and incase other databases are used in the future.

    Args:
        name_of_nuclide (str): The name of the nuclide in the form PU243.
        df_containing_nuclides (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.
        nuclide_in_ORION (bool): Boolean operator that is defaulted to False. Only set to True if the nuclide exits in the ORION database.

    Returns:
        orion_id_of_nuclide (int or str): The ORION ID of the nuclide. This is the index of the array containing the specified nuclide, unless the nuclide does not exist in the ORION database, then it is set to the given string.
        nuclide_in_ORION (bool): Boolean operator that is defaulted to False. Only set to True if the nuclide exits in the ORION database.
    """
    # The nuclide names are give as 12-MG33 in the ORION database so this removes the N- prefix
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in df_containing_nuclides]

    # Iterates through each row of the array containing the ORION nuclides
    for row in df_nuclide_without_prefix:
        # If the current row matches the name of the examined nuclide then the nuclide exists in the ORION database and the loop is exited
        if row == name_of_nuclide:
            orion_id_of_nuclide = df_nuclide_without_prefix.index(row) # ORION ID is set to the index of the row containing the examined nuclide
            nuclide_in_ORION = True
            break
    
    # If the boolean operator remains false then the nuclide does not exist in the ORION database
    if nuclide_in_ORION == False:
        orion_id_of_nuclide = 'Non Existant'
    
    return orion_id_of_nuclide, nuclide_in_ORION

def if_parent_daughter_exist_in_orion(zaid_of_parent_daughter, df_containing_nuclides, parent_daughter_exists_in_orion):
    """
    Determines if the parent/daughter nuclide exist in the ORION database.
    This is done because the U241 nuclide exists in the ENDF database but not the ORION one. Hence, this makes sure any nuclides with a parent/daughter of U241 do not have the ZAID of U241 (922410) as an output.
    This is the only case of this but the process has been generalised for validation and incase other databases are used in the future.

    Args:
        zaid_of_parent_daughter (int): The ZAID of the possible parent/daughter nuclide, in the form ZZAAAI (I=0).
        df_containing_nuclides (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.
        parent_daughter_exists_in_orion (bool): Boolean operator that is defaulted to False. Only set to True if the parent/daughter nuclide exits in the ORION database.

    Returns:
        parent_daughter_exists_in_orion (bool): Boolean operator that is defaulted to False. Only set to True if the parent/daughter nuclide exits in the ORION database.
    """
    # Gather the parent/daughter nuclide information (name, Z, A)
    nuclide_details = build_nuclide_name(zaid_of_parent_daughter)
    parent_daughter_name = nuclide_details[0]
    
    # The nuclide names are give as 12-MG33 in the ORION database so this removes the N- prefix
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in df_containing_nuclides]
    parent_daughter_exists_in_orion = False
    # Iterates through each row of the array containing the ORION nuclides
    for row in df_nuclide_without_prefix:
        # If the current row matches the name of the parent/daughter nuclide then the parent/daughter nuclide exists in the ORION database and the loop is exited
        if parent_daughter_name == row:
            parent_daughter_exists_in_orion = True
            break
    
    return parent_daughter_exists_in_orion

def metastable_parent(zaid_metastable_precursor, nuclide_in_dataframe):
    """_summary_

    Args:
        zaid_metastable_precursor (_type_): _description_
        nuclide_in_dataframe (_type_): _description_

    Returns:
        _type_: _description_
    """
    metastable_list_dir = f'{args.reaction_data}/NNL_metastable_precursors_and_parents.csv'
    metastable_df = pd.read_csv(metastable_list_dir, header = 0)
    
    # Filters through the 'ZAID' column of the pandas dataframe containing the metastable nuclides and extracts the row that has the same ZAID as the precursor nuclide
    examined_precursor = metastable_df[(metastable_df['ZAID'] == zaid_metastable_precursor)]
    
    # Creates a pandas series(array) that contains the parent ZAIDs of the metastable precursor
    metastable_parent_list = pd.Series(examined_precursor.iloc[:, 3: 15].values.flatten().tolist())
    # The parents with 'NA' in the .csv are saved as nan in the array, so these are converted back to 'NA' by lambda
    # All of the ZAIDs are converted from floats to integers
    metastable_parent_list = metastable_parent_list.apply(lambda x: 'NA' if pd.isna(x) else int(x))
    # Loops through each row in the array of metastable parents and checks if they exist in the ORION database
    # Enumerates the array so that the index of the parent is found, so that the relevant parent can be inspected
    for i, row in enumerate(metastable_parent_list):
        metastable_parent_exist = False
        # The ORION database check is skipped out if the nuclide has parent ZAID == 'NA'
        if row != 'NA':
            metastable_parent_exist = if_parent_daughter_exist_in_orion(row, nuclide_in_dataframe, metastable_parent_exist)
            if metastable_parent_exist == False:
                metastable_parent_list[i] = 'NA'

    return metastable_parent_list

def gs_pecursor_es_parent(empty_parent_list, zaid_gs_precursor, num_parents):
    """_summary_

    Args:
        empty_parent_list (_type_): _description_
        zaid_gs_precursor (_type_): _description_
        num_parents (_type_): _description_

    Returns:
        _type_: _description_
    """
    reactions_of_es_precursor_dir = f'{args.reaction_data}/NNL_metastable_precursors_and_daughters.csv'
    df_es_precursors = pd.read_csv(reactions_of_es_precursor_dir, header = 0)

    # Checks to see if the ZAID of a ground state nuclide is in the list of daughter nuclides from metastable precursors
    examined_daughter_16 = df_es_precursors[(df_es_precursors['Daughter 16 ZAID'] == zaid_gs_precursor)]
    # If the ground state ZAID of the daughter exists in the dataframe then the metastable precursor ZAID is assigned to the relevant index in the list of parents
    if not examined_daughter_16.empty:
        es_precursor = examined_daughter_16['Precursor ZAID'].iloc[0]
        if (es_precursor % 10 == 1):
            empty_parent_list[4] = es_precursor # First excited state of parent
            num_parents += 1
        elif (es_precursor % 10 == 2):
            empty_parent_list[5] = es_precursor # Second excited state of parent
            num_parents += 1

    examined_daughter_17 = df_es_precursors[(df_es_precursors['Daughter 17 ZAID'] == zaid_gs_precursor)]
    if not examined_daughter_17.empty:
        es_precursor = examined_daughter_17['Precursor ZAID'].iloc[0]
        if (es_precursor % 10 == 1):
            empty_parent_list[6] = es_precursor
            num_parents += 1
        elif (es_precursor % 10 == 2):
            empty_parent_list[7] = es_precursor
            num_parents += 1

    examined_daughter_102 = df_es_precursors[(df_es_precursors['Daughter 102 ZAID'] == zaid_gs_precursor)]
    if not examined_daughter_102.empty:
        es_precursor = examined_daughter_102['Precursor ZAID'].iloc[0]
        if (es_precursor % 10 == 1):
            empty_parent_list[8] = es_precursor
            num_parents += 1
        elif (es_precursor % 10 == 2):
            empty_parent_list[9] = es_precursor
            num_parents += 1
    
    examined_daughter_103 = df_es_precursors[(df_es_precursors['Daughter 103 ZAID'] == zaid_gs_precursor)]
    if not examined_daughter_103.empty:
        es_precursor = examined_daughter_103['Precursor ZAID'].iloc[0]
        if (es_precursor % 10 == 1):
            empty_parent_list[10] = es_precursor
            num_parents += 1
        elif (es_precursor % 10 == 2):
            empty_parent_list[11] = es_precursor
            num_parents += 1
    
    return empty_parent_list, num_parents

def is_isomer_in_NNL(metastable_zaid, exists_in_NNL):
    """_summary_

    Args:
        metastable_zaid (_type_): _description_
        exists_in_NNL (_type_): _description_

    Returns:
        _type_: _description_
    """
    metastable_list_dir = f'{args.reaction_data}/NNL_metastable_precursors_and_parents.csv'
    metastable_df = pd.read_csv(metastable_list_dir, header = 0)
    
    # Filters through the 'ZAID' column of the pandas dataframe containing the metastable nuclides and extracts the row that has the same ZAID as the precursor nuclide
    examined_precursor = metastable_df[(metastable_df['ZAID'] == metastable_zaid)]

    if examined_precursor.empty:
        exists_in_NNL = False

    return exists_in_NNL
# -------------  Main Code  ------------ #
args = parse_arguments()
# Define the Element Symbols
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md']

# Opens the csv file containing the ENDF and JEFF reaction data for each nuclide and stores it in an array
reaction_file = f'{args.cross_section_data}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_all_reactiondata.csv'
reaction_data = np.genfromtxt(reaction_file, comments = '%', delimiter = ',', skip_header = 1)

# Opens the excel file containing the ORION IDs for each nuclide and loads it in a pandas dataframe
ORION_ID_dir = f'{args.reaction_data}/orion_nuclides_list.xlsx'
df_ORION_ID = pd.read_excel(ORION_ID_dir, header = None)
# Columns in excel file do not contain headers so create them here
df_ORION_ID.columns = ['Nuclide Name', 'Buffer Mass']
nuclide_in_df = df_ORION_ID[df_ORION_ID.columns[0]]

nuclide_data =[]
previous_nuclide = None
grouped_endf_jeff_dictionary = {}
# Iterates through each line of the ENDF and JEFF data
for row in reaction_data:
    # Sets the nuclide being examined to the current row in the loop
    current_nuclide = row[0]
    if current_nuclide == previous_nuclide:
        # If the ZAID in the examined row of the .csv file is the same as the previous row then the data in the current row is appended to the same group in the dictionary as the previous row
        grouped_endf_jeff_dictionary[current_nuclide].append(row)
    else:
        if current_nuclide not in grouped_endf_jeff_dictionary:
            # If the ZAID in the examined row is not a group in the dictionary then a group is created for this nuclide
            grouped_endf_jeff_dictionary[current_nuclide] = []
        # The current row is appended to the group for the new nuclide
        grouped_endf_jeff_dictionary[current_nuclide].append(row)

    previous_nuclide = current_nuclide


# Iterates through each dictionary group (key, which is the ZAID of the nuclide) and the reaction data for each reaction
for key, reactions in grouped_endf_jeff_dictionary.items():
    nuclide_characteristics = build_nuclide_name(key)
    nuclide_name = nuclide_characteristics[0]
    
    # Check whether the nuclide is in the ORION database at all
    orion_id, exists_in_ORION = find_ORION_ID(nuclide_name, nuclide_in_df, nuclide_in_ORION = False)
    if exists_in_ORION == True:

        is_metastable = key % 10
        if is_metastable != 0:

            isomer_exists_in_NNL = True
            isomer_exists_in_NNL = is_isomer_in_NNL(key, isomer_exists_in_NNL)
            # Check to see if the isomer exists in the list of metastable isomers in the NNL database, eg Ag117(1) is in JEFF but not NNL
            if isomer_exists_in_NNL == True:
                parent_list = metastable_parent(key, nuclide_in_df)
                parent_16 = parent_list[0]
                parent_17 = parent_list[1]
                parent_102 = parent_list[2]
                parent_103 = parent_list[3]
                NumberParents = 0
                for parent in parent_list:
                    if parent != 'NA':
                        NumberParents = NumberParents + 1
            else:
                # If the isomer is not in the NNL database then it's parents are all set to 'NA' and NumberParents = 0
                NumberParents = 0
                parent_list = pd.Series(['NA'] * 12)
                # Resets the values of these variables. This is necessary as if the precursor nuclide prior has parents in any of these reactions then a metastable nuclide that exists in JEFF but not NNL must have all parents assigned to 'NA'
                parent_16 = parent_17 = parent_102 = parent_103 = 'NA'

        else:
            NumberParents, parent_16, parent_17, parent_102, parent_103 = calculate_parent(nuclide_characteristics[1],
                                                                                            nuclide_characteristics[2],
                                                                                            nuclide_in_df)
            parent_list = pd.Series(['NA'] * 12)
            parent_list, NumberParents = gs_pecursor_es_parent(parent_list, key, NumberParents)

        reaction_16_exists = reaction_17_exists = reaction_18_exists = reaction_102_exists = reaction_103_exists = False
        # Iterates through each column in the arrays containing the reaction data
        for column in reactions:
            # If the MT value of the reaction matches the if statement then the daughter and parent IDs are extracted
            if column[1] == 16:
                reaction_16_exists = True
                daughter_16 = int(column[2]) # Change ID from type numpy.float64 to typeinteger
                daughter_exists = find_daughter(daughter_16, nuclide_in_df)
                if daughter_exists == False:
                    daughter_16 = 'NA'

            if column[1] == 17:
                reaction_17_exists = True
                daughter_17 = int(column[2])
                daughter_exists = find_daughter(daughter_17, nuclide_in_df)
                if daughter_exists == False:
                    daughter_17 = 'NA'

            if column[1] == 18:
                reaction_18_exists = True
                daughter_18 = int(column[2])
                daughter_exists = find_daughter(daughter_18, nuclide_in_df)
                if daughter_exists == False:
                    daughter_18 = 'NA'

            if column[1] == 102:
                reaction_102_exists = True
                daughter_102 = int(column[2])
                daughter_exists = find_daughter(daughter_102, nuclide_in_df)
                if daughter_exists == False:
                    daughter_102 = 'NA'

            if column[1] == 103:
                reaction_103_exists = True
                daughter_103 = int(column[2])
                daughter_exists = find_daughter(daughter_103, nuclide_in_df)
                if daughter_exists == False:
                    daughter_103 = 'NA'

        # If the reaction does not exist in the JEFF database then the daughter ZAID is set to 'NA'
        if reaction_16_exists == False:
            daughter_16 = 'NA'
        if reaction_17_exists == False:
            daughter_17 = 'NA'
        if reaction_18_exists == False:
            daughter_18 = 'NA'
        if reaction_102_exists == False:
            daughter_102 = 'NA'
        if reaction_103_exists == False:
            daughter_103 = 'NA'

        # Calculate the number of reactions for which the nuclide is the parent
        daughter_ids_array = [daughter_16, daughter_17, daughter_18, daughter_102, daughter_103]
        NumberReactions = 0
        for value in daughter_ids_array:
            if value != 'NA':
                NumberReactions = NumberReactions + 1
        if key == 471181:
            print('Hi')
            print(parent_list)
            print(parent_102)
        nuclide_data.append({'Nuclide Name': nuclide_name, 
                                    'Nuclide ZAID': int(key),
                                    'ORION ID': orion_id, 
                                    'NumberReactions': NumberReactions, 
                                    'Reaction 16 Daughter ZAID': daughter_16, 
                                    'Reaction 17 Daughter ZAID': daughter_17, 
                                    'Reaction 18 Daughter ZAID': daughter_18, 
                                    'Reaction 102 Daughter ZAID': daughter_102, 
                                    'Reaction 103 Daughter ZAID': daughter_103, 
                                    'Number Parents': NumberParents, 
                                    'Reaction 16 Parent ZAID': parent_16, 
                                    'Reaction 17 Parent ZAID': parent_17, 
                                    'Reaction 102 Parent ZAID': parent_102, 
                                    'Reaction 103 Parent ZAID': parent_103,
                                    'Reaction 16 Parent(M) ZAID': parent_list[4],
                                    'Reaction 16 Parent(2M) ZAID': parent_list[5],
                                    'Reaction 17 Parent(M) ZAID': parent_list[6],
                                    'Reaction 17 Parent(2M) ZAID': parent_list[7],
                                    'Reaction 102 Parent(M) ZAID': parent_list[8],
                                    'Reaction 102 Parent(2M) ZAID': parent_list[9],
                                    'Reaction 103 Parent(M) ZAID': parent_list[10],
                                    'Reaction 103 Parent(2M) ZAID': parent_list[11]})

        print('Appending data for', nuclide_name)





# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(nuclide_data)

# Save the whole array to a csv file
file_name = f'{args.path}/ZAID_results.csv'
df.to_csv(file_name, index = False)

