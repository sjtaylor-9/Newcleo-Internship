"""
Finding_parent_daughter_nuclides.py

This code is responsible for finding the ZAIDs of the parent and daughter nuclides, which are needed for the MPR file. The ZAIDs of the parent and daughter nuclides for each reaction are outputted to a .csv file.
In addition, the ORION IDs NumberParents, and NumberReactions are also written to the .csv file. These are also all neded for the MPR file.
The code calculates the parent/daughter ZAIDs for each nuclide and then determines if these exist in both the ENDF and ORION databases. All of the reaction information for the nuclide is then appended to a pandas dataframe and written to a .csv file once all nuclides that exist in both the ENDF and ORION databases have been examined.

Author: Sam Taylor (sam.taylor@newcleo.com)
Last Edited: 22/11/2024
"""
# ---------  Import Libraries  --------- #
import numpy as np
import endf
import os
import pandas as pd
# -------------  Functions  ------------ #
def build_nuclide_name(nuclide_string):
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
    # Set Z to the first digit if the second digit is '0' and the nuclide ID number is less than 9999
    if nuclide < 10000:
        atomic_number = int(nuclide_string[0])
    else:
        # Otherwise, set Z to the first two digits
        atomic_number = int(nuclide_string[:2])

    # Mass number is set to figure left after removing first 2 digits
    mass_number = nuclide_string[2:]

    # Remove leading zeros from A; if A is empty after stripping, set it to 0
    mass_number = int(mass_number.lstrip('0')) if mass_number.lstrip('0') else 0

    # Nuclide name in form PU243
    nuclear_notation = element_symbols[atomic_number-1].upper() + str(mass_number)
    
    return nuclear_notation, atomic_number, mass_number

def build_nuclide_ID(atomic_number, mass_number):
    """
    Constructs the ZAID of the nuclide from it's atomic and mass numbers. The ZAIDs are used to identify the nuclides.
    The ZAID is in the form ZZAAI, where Z is the atomic number and A is the mass number. I > 0 identifies that an isotope is metastable. In this code I is set to 0 since I > 0 nuclides are not in the ENDF databse.

    Args:
        atomic_number (int): The atomic number of the nuclide, Z.
        mass_number (int): The mass number of the nuclide, A.

    Returns:
        ID (int): The ZAID number of the nuclide in the form ZZAAAI (I is set to 0, I > 0 for metastable isotopes).
    """
    # ZAID calculated as 1000*Z + A
    if mass_number < 10:
        ID = str(atomic_number) + "00" + str(mass_number)
    else:
        if mass_number < 100:
            ID = str(atomic_number) + "0" + str(mass_number)
        else:
            ID = str(atomic_number) + str(mass_number)
    
    # Nuclear IDs in form ZZAAAI, where I = 0,1,2. I > 0 is for metastable istotopes. For now ignoring these as not in ENDF database.
    ID = str(ID) + "0"
    
    return int(ID)

def calculate_daughter(atomic_number, mass_number, nuclide_in_dataframe, reaction):
    """
    Determines if the daughter nuclide of a given reaction exists in both the ENDF and ORION databases and if so then the ZAID of the daughter nuclide is calculated.
    If the daughter nuclide does not exist in either databse then the ZAID is set to 'NA'.
    
    MT IDs of reactions under consideration are 16, 17, 18, 102, 103.
    Reaction 16 is (n,2n) and so one initial neutron absorbed causing 2 neutrons to be emitted, so A decreases by 1.
    Reaction 17 is (n,3n) and so one initial neutron absorbed causing 3 neutrons to be emitted, so A decreases by 2.
    Reaction 18 is (n,fission) and so one initial neutron causes fission to generate fission products.
    Reaction 102 is (n,gamma) and so a neutron is absorbed causing nuclear excitation, which results in gamma emission, so A increases by 1.
    Reaction 103 is (n,p) and so a neutron is absorbed causing proton emission, so Z decreases by 1 and A is constant.

    Args:
        atomic_number (int): The atomic number (Z) of the nuclide under examination (i.e. the parent of the possible daughter nuclide).
        mass_number (int): The mass number (A) of the nuclide under examination (i.e. the parent of the possible daughter nuclide).
        nuclide_in_dataframe (pandas dataframe): The pandas dataframe storing the list of nuclides in the ORION database in the order of ascending ORION IDs.
        reaction (int): The MT value of the given reaction (16, 17, 18, 102, 103).

    Returns:
        daughter_ID_MT (int or str): Where (MT = 16, 17, 18, 102, 103). If the daughter nuclide of the specified reaction exists in both the ENDF and ORION databases, then the ZAID is outputted. However, if it does not then the ZAID is set to 'NA'. 
                                     If the reaction has MT = 18, then the examined nuclide has fissioned and so the ZAID of the daughter is set to 0, as too difficult to determine the fission product.
    """
    
    if reaction == 16:
        # If nuclide is parent of reaction 16 (n,2n) then daughter should have A-1
        daughter_A_16 = mass_number - 1
        daughter_ID_16 = build_nuclide_ID(atomic_number, daughter_A_16)
        # Checks to see if the daughter nuclide exists in the database
        daughter_exists = find_daughter(daughter_ID_16, nuclide_in_dataframe)
        if daughter_exists == False:
            daughter_ID_16 = 'NA'
        
        return daughter_ID_16
    
    elif reaction == 17:
        # If nuclide is parent of reaction 17 (n,3n) then daughter should have A-2
        daughter_A_17 = mass_number - 2
        daughter_ID_17 = build_nuclide_ID(atomic_number, daughter_A_17)
        daughter_exists = find_daughter(daughter_ID_17, nuclide_in_dataframe)
        if daughter_exists == False:
            daughter_ID_17 = 'NA'
        
        return daughter_ID_17
    
    elif reaction == 18:
        # If nuclide is parent of reaction 18 (n,fission) then daughter ID is 0
        daughter_ID_18 = 0
        return daughter_ID_18
    
    elif reaction == 102:
        # If nuclide is parent of reaction 102 then daughter should have A+1
        daughter_A_102 = mass_number + 1
        daughter_ID_102 = build_nuclide_ID(atomic_number, daughter_A_102)
        daughter_exists = find_daughter(daughter_ID_102, nuclide_in_dataframe)
        if daughter_exists == False:
            daughter_ID_102 = 'NA'
        
        return daughter_ID_102
    
    elif reaction == 103:
        # If nuclide is parent of reaction 103 then daughter should have same A and Z-1
        daughter_Z_103 = atomic_number - 1
        daughter_ID_103 = build_nuclide_ID(daughter_Z_103, mass_number)
        daughter_exists = find_daughter(daughter_ID_103, nuclide_in_dataframe)
        if daughter_exists == False:
            daughter_ID_103 = 'NA'
        
        return daughter_ID_103

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
        if row[0] == parent_id and row[3] == MT_value:
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
        if row[0] == daughter_id:
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
    nuclide_details = build_nuclide_name(str(zaid_of_parent_daughter))
    parent_daughter_name = nuclide_details[0]
    # Removes the 0 at the end of the parent/daughter nuclide name, which is due to the input being in the form ZZAAAI rather than ZZAAA
    parent_daughter_name = parent_daughter_name[:-1]
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
# -------------  Main Code  ------------ #

# Define the Element Symbols. Eisteinium included to prevent list index issues when calculating the parent nuclide name of Cf.
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']

# Path to Lib80x files, which are the cross-sections of reactions in the ENDF database
lib_dir = r'/Users/sam/Documents/NewcleoInternship/Lib80x/Lib80x/'

# Opens the csv file containing the reaction data for each nuclide and stores it in an array
reaction_data_dir = r'/Users/sam/Documents/NewcleoInternship/'
reactor_type = r'LFR30_MPR_reactiondata.csv'
reaction_file = reaction_data_dir + reactor_type
reaction_data = np.genfromtxt(reaction_file, comments = '%', delimiter = ',')

# Opens the excel file containing the ORION IDs for each nuclide and loads it in a pandas dataframe
ORION_ID_dir = r'/Users/sam/Documents/NewcleoInternship/orion_nuclides_list.xlsx'
df_ORION_ID = pd.read_excel(ORION_ID_dir, header = None)
# Columns in excel file do not contain headers so create them here
df_ORION_ID.columns = ['Nuclide Name', 'Buffer Mass']
nuclide_in_df = df_ORION_ID[df_ORION_ID.columns[0]]

nuclide_data =[]
nuclide = 1001
while nuclide < 98255: # Last nuclide is 254Cf
    nuclide_str = str(nuclide)

    nuclide_name, Z, A = build_nuclide_name(nuclide_str)
    nuclide_ID = build_nuclide_ID(Z, A)
    # Build infile name based on chosen nuclide
    infile_str = (lib_dir + element_symbols[Z - 1] + "/" + nuclide_str + ".802nc")

    # Check whether the nuclide is in the database at all by validating if ENDF file exists
    if os.path.isfile(infile_str):
        
        # Check whether the nuclide is in the ORION database at all
        orion_id, exists_in_ORION = find_ORION_ID(nuclide_name, nuclide_in_df, nuclide_in_ORION = False)
        # The nuclide data is only appended to the array if the nuclide exists in both ENDF and ORION. For eg: U241 exists in the ENDF database but not in the ORION database
        if exists_in_ORION == True:
            table = endf.ace.get_table(infile_str)

            # Check whether the database has the required reactions for the nuclide and if so calculate the daughter nuclide ID in the form ZZAAAI
            # If the reaction does not exist then the ID is set to NA
            if 16 in (table.interpret().reactions).keys():
                daughter_16 = calculate_daughter(Z, A, nuclide_in_df, reaction = 16)
            else:
                daughter_16 = 'NA'

            if 17 in (table.interpret().reactions).keys():
                daughter_17 = calculate_daughter(Z, A, nuclide_in_df, reaction = 17)
            else:
                daughter_17 = 'NA'

            if 18 in (table.interpret().reactions).keys():
                daughter_18 = calculate_daughter(Z, A, nuclide_in_df, reaction = 18)
            else:
                daughter_18 = 'NA'

            if 102 in (table.interpret().reactions).keys():
                daughter_102 = calculate_daughter(Z, A, nuclide_in_df, reaction = 102)
            else:
                daughter_102 = 'NA'

            if 103 in (table.interpret().reactions).keys():
                daughter_103 = calculate_daughter(Z, A, nuclide_in_df, reaction = 103)
            else:
                daughter_103 = 'NA'

            # Calculate the number of reactions for which the nuclide is the parent
            daughter_ids_array = [daughter_16, daughter_17, daughter_18, daughter_102, daughter_103]
            NumberReactions = 0
            for value in daughter_ids_array:
                if value != 'NA':
                    NumberReactions = NumberReactions + 1

            NumberParents, parent_16, parent_17, parent_102, parent_103 = calculate_parent(Z, A, nuclide_in_df)
            
            nuclide_data.append({'Nuclide Name': nuclide_name, 
                                    'Nuclide ZAID': nuclide_ID,
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
                                    'Reaction 103 Parent ZAID': parent_103,})

            print('Appending data for', nuclide_name)

    nuclide = nuclide + 1

# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(nuclide_data)

# Save the whole array to a csv file
file_name = 'ZAID_results.csv'
df.to_csv(file_name, index = False)

