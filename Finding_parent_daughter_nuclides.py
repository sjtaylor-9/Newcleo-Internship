import numpy as np
import endf
import os
import pandas as pd

def build_nuclide_name(nuclide_string):
    
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
    """_summary_

    Args:
        atomic_number (_type_): _description_
        mass_number (_type_): _description_

    Returns:
        _type_: _description_
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
    MT IDs of reactions under consideration are 16, 17, 18, 102, 103.
    Reaction 16 is (n,2n) and so one initial neutron absorbed causing 2 neutrons to be emitted, so A decreases by 1.
    Reaction 17 is (n,3n) and so one initial neutron absorbed causing 3 neutrons to be emitted, so A decreases by 2.
    Reaction 18 is (n,fission) and so one initial neutron causes fission to generate fission products.
    Reaction 102 is (n,gamma) and so a neutron is absorbed causing nuclear excitation, which results in gamma emission, so A increases by 1.
    Reaction 103 is (n,p) and so a neutron is absorbed causing proton emission, so Z decreases by 1 and A is constant.

    Args:
        atomic_number (_type_): _description_
        mass_number (_type_): _description_
        reaction (_type_): _description_

    Returns:
        _type_: _description_
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
    """_summary_

    Args:
        atomic_number (_type_): _description_
        mass_number (_type_): _description_

    Returns:
        _type_: _description_
    """
    # If nuclide is daughter of reaction 16 (n,2n) then parent should have A+1
    parent_A_16 = mass_number + 1
    parent_16 = build_nuclide_ID(atomic_number, parent_A_16)
    
    # If nuclide is daughter of reaction 17 (n,3n) then parent should have A+2
    parent_A_17 = mass_number + 2
    parent_17 = build_nuclide_ID(atomic_number, parent_A_17)

    # If nuclide is daughter of reaction 18 (n,fission) then parent ID is 0

    # If nuclide is daughter of reaction 102 then parent should have A-1
    parent_A_102 = mass_number - 1
    parent_102 = build_nuclide_ID(atomic_number, parent_A_102)

    # If nuclide is daughter of reaction 103 then parent should have same A and Z+1
    parent_Z_103 = atomic_number + 1
    parent_103 = build_nuclide_ID(parent_Z_103, mass_number)

    # Checks to see if the parent nuclide is in the database by comparing the possible parent ID with the ones in the first column of the csv file
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
    """_summary_

    Args:
        NumberParents (_type_): _description_
        parent_id (_type_): _description_
        MT_value (_type_): _description_

    Returns:
        _type_: _description_
    """
    parent_exists = False
    for row in reaction_data:
        if row[0] == parent_id and row[3] == MT_value:
            parent_exists = True
            parent_exists = if_parent_daughter_exist_in_orion(parent_id, df_containing_nuclides, parent_exists)
            if parent_exists == True:
                NumberParents = NumberParents + 1
                break
            
    return NumberParents, parent_exists

def find_daughter(daughter_id, df_containing_nuclides):
    """_summary_

    Args:
        daughter_id (_type_): _description_

    Returns:
        _type_: _description_
    """
    daughter_exists = False

    for row in reaction_data:
        if row[0] == daughter_id:
            daughter_exists = True
            daughter_exists = if_parent_daughter_exist_in_orion(daughter_id, df_containing_nuclides, daughter_exists)
            break
    
    return daughter_exists

def find_ORION_ID(name_of_nuclide, df_containing_nuclides, nuclide_in_ORION):
    """_summary_

    Args:
        name_of_nuclide (_type_): _description_
        df_containing_nuclides (_type_): _description_
        nuclide_in_ORION (_type_): _description_

    Returns:
        _type_: _description_
    """
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in df_containing_nuclides]

    for row in df_nuclide_without_prefix:
        if row == name_of_nuclide:
            orion_id_of_nuclide = df_nuclide_without_prefix.index(row)
            nuclide_in_ORION = True
            break
    
    if nuclide_in_ORION == False:
        orion_id_of_nuclide = 'Non Existant'
    
    return orion_id_of_nuclide, nuclide_in_ORION

def if_parent_daughter_exist_in_orion(zaid_of_parent_daughter, df_containing_nuclides, parent_daughter_exists_in_orion):
    """_summary_

    Args:
        zaid_of_parent_daughter (_type_): _description_
        df_containing_nuclides (_type_): _description_
        parent_daughter_exists_in_orion (_type_): _description_

    Returns:
        _type_: _description_
    """
    nuclide_details = build_nuclide_name(str(zaid_of_parent_daughter))
    parent_daughter_name = nuclide_details[0]
    parent_daughter_name = parent_daughter_name[:-1]
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in df_containing_nuclides]
    parent_daughter_exists_in_orion = False
    for row in df_nuclide_without_prefix:
        if parent_daughter_name == row:
            parent_daughter_exists_in_orion = True
            break
    
    return parent_daughter_exists_in_orion

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

