# Imports the necessary Python packages
import numpy as np
import pandas as pd

def calculate_parent_orion(parent_zaid, orion_data):
    """_summary_

    Args:
        orion (_type_): _description_
        reaction (_type_): _description_

    Returns:
        _type_: _description_
    """
    parent_name = build_parent_daughter_name(parent_zaid)
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in orion_data]

    for row in df_nuclide_without_prefix:
        if row == parent_name[0]:
            parent_orion = df_nuclide_without_prefix.index(row)
            break
    return parent_orion

def calculate_daughter_orion(daughter_zaid, orion_data):
    """_summary_

    Args:
        daughter_zaid (_type_): _description_
        orion_data (_type_): _description_

    Returns:
        _type_: _description_
    """
    daughter_name = build_parent_daughter_name(daughter_zaid)
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in orion_data]

    for row in df_nuclide_without_prefix:
        if row == daughter_name[0]:
            daughter_orion = df_nuclide_without_prefix.index(row)
            break

    return daughter_orion

def build_parent_daughter_name(id_number):
    """_summary_

    Args:
        id_number (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Format the ZAIDs correctly, in both str and int types, such that .0 is removed and the trailing (I=0) is removed as well
    id_number = int(id_number)
    id_string = str(id_number)
    if id_string.endswith('0'):
        id_string = id_string[:-1]
    id_number = int(id_string)

    # Set Z to the first digit if the second digit is '0' and the nuclide ID number is less than 9999
    if id_number < 10000:
        atomic_number = int(id_string[0])
    else:
        # Otherwise, set Z to the first two digits
        atomic_number = int(id_string[:2])
   
    # Mass number is set to figure left after removing first 2 digits
    mass_number = id_string[2:]

    # Remove leading zeros from A; if A is empty after stripping, set it to 0
    mass_number = int(float(mass_number.lstrip('0'))) if mass_number.lstrip('0') else 0

    # Nuclide name in form PU243
    nuclear_notation = element_symbols[atomic_number-1].upper() + str(mass_number)
    return nuclear_notation, atomic_number, mass_number

def get_cross_section(zaid_in_MPR, reaction):
    """_summary_

    Args:
        zaid_in_MPR (_type_): _description_
        reaction (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Opens the csv file containing the cross-section data for each nuclide and stores it in an array
    #xs_data_dir = r'C:\Users\sam.taylor\OneDrive - Newcleo\Documents\Modelling_LFR\Generating_MPR_file'
    xs_data_dir = r'/Users/sam/Documents/NewcleoInternship/'
    reactor_type = r'LFR30_MPR_reactiondata.csv'
    xs_file = xs_data_dir + reactor_type
    xs_data = np.genfromtxt(xs_file, comments = '%', delimiter = ',')
    xs = 0
    # Iterates through the csv file until the required reaction is found and the cross-section is outputted
    for row in xs_data:
        if row[0] == zaid_in_MPR and row[3] == reaction:
            xs = row[4]
            break
        
    return xs

def remove_trailing_zeroes(nuclide_id_with_zeroes):
    nuclide_id_without_zeroes = str(nuclide_id_with_zeroes).rstrip('.0')
    nuclide_id_without_zeroes = nuclide_id_without_zeroes + "0"
    return nuclide_id_without_zeroes

def save_to_MPR(ORION, nuclide_ID, name, NumberReactions, daughter_16, daughter_17, daughter_18, daughter_102, daughter_103, NumberParents, parent_16, parent_17, parent_102, parent_103, orion_dataframe):
    """_summary_

    Args:
        ORION (_type_): _description_
        nuclide_ID (_type_): _description_
        name (_type_): _description_
        NumberReactions (_type_): _description_
        daughter_16 (_type_): _description_
        daughter_17 (_type_): _description_
        daughter_18 (_type_): _description_
        daughter_102 (_type_): _description_
        daughter_103 (_type_): _description_
        NumberParents (_type_): _description_
        parent_16 (_type_): _description_
        parent_17 (_type_): _description_
        parent_102 (_type_): _description_
        parent_103 (_type_): _description_
    """
    # The ORION and ZAIDs are right aligned with total widths of 5 and 8, respectively, and the nuclide name is left aligned with a total width of 11
    width_orion = 5
    aligned_orion = "{:>{width}}".format(ORION, width = width_orion)
    width_nuclide = 8
    aligned_nuclide = "{:>{width}}".format(nuclide_ID, width = width_nuclide)
    width_name = 11
    aligned_name = "{:<{width}}".format(name, width = width_name)
    
    # NumberReactions and NumberParents should always be < 10 but have right aligned as validation
    width_num_reactions = 5
    aligned_reactions = "{:>{width}}".format(NumberReactions, width = width_num_reactions)
    width_num_parents = 4
    aligned_parents = "{:>{width}}".format(NumberParents, width = width_num_parents)
    file.write(f"{aligned_orion}{aligned_nuclide}  {aligned_name}NumberReactions{aligned_reactions} NumberParents{aligned_parents}\n")
    
    # Adds the parent IDs in the required format
    # Order of reactions in increasing ORION IDs is 102, 16, 17, 103, no ORION ID for reaction 18
    # Text for both IDs must be right aligned both with total width of 12
    width_parent_orion = 12
    width_parent_zaid = 12
    if not np.isnan(parent_102):
        parent_orion_102 = calculate_parent_orion(parent_102, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_102, width = width_parent_orion)
        parent_102 = remove_trailing_zeroes(parent_102)
        aligned_zaid = "{:>{width}}".format(parent_102, width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_16):
        parent_orion_16 = calculate_parent_orion(parent_16, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_16, width = width_parent_orion)
        parent_16 = remove_trailing_zeroes(parent_16)
        aligned_zaid = "{:>{width}}".format(parent_16, width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_17):
        parent_orion_17 = calculate_parent_orion(parent_17, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_17, width = width_parent_orion)
        parent_17 = remove_trailing_zeroes(parent_17)
        aligned_zaid = "{:>{width}}".format(parent_17, width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_103):
        parent_orion_103 = calculate_parent_orion(parent_103, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_103, width = width_parent_orion)
        parent_103 = remove_trailing_zeroes(parent_103)
        aligned_zaid = "{:>{width}}".format(parent_103, width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    
    # Adds the relevant data for the daughter nuclides in the required format
    # The ORION ID of the daughter nuclide must be right aligned with a total width of 7 (6 preceeding spaces for 1 digit) and ZAID must be right aligned with a total width of 8
    width_orion = 7
    width_zaid = 8
    if not np.isnan(daughter_16):
        daughter_orion_16 = calculate_daughter_orion(daughter_16, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_16, width = width_orion)
        daughter_16 = remove_trailing_zeroes(daughter_16)
        aligned_zaid = "{:>{width}}".format(daughter_16, width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  16 (n,2n)\n")
        file.write(f"{get_cross_section(nuclide_ID, reaction = 16)}\n")
    
    if not np.isnan(daughter_17):
        daughter_orion_17 = calculate_daughter_orion(daughter_17, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_17, width = width_orion)
        daughter_17 = remove_trailing_zeroes(daughter_17)
        aligned_zaid = "{:>{width}}".format(daughter_17, width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  17 (n,3n)\n")
        file.write(f"{get_cross_section(nuclide_ID, reaction = 17)}\n")
    
    if not np.isnan(daughter_18):
        daughter_orion_18 = -1 # ORION ID for fission product is -1
        aligned_orion = "{:>{width}}".format(daughter_orion_18, width = width_orion)
        daughter_18 = remove_trailing_zeroes(daughter_18)
        aligned_zaid = "{:>{width}}".format(daughter_18, width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  18 (n,fission)\n")
        file.write(f"{get_cross_section(nuclide_ID, reaction = 18)}\n")

    if not np.isnan(daughter_102):
        daughter_orion_102 = calculate_daughter_orion(daughter_102, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_102, width = width_orion)
        daughter_102 = remove_trailing_zeroes(daughter_102)
        aligned_zaid = "{:>{width}}".format(daughter_102, width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  102 (n,gamma)\n")
        file.write(f"{get_cross_section(nuclide_ID, reaction = 102)}\n")

    if not np.isnan(daughter_103):
        daughter_orion_103 = calculate_daughter_orion(daughter_103, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_103, width = width_orion)
        daughter_103 = remove_trailing_zeroes(daughter_103)
        aligned_zaid = "{:>{width}}".format(daughter_103, width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  103 (n,proton)\n")
        file.write(f"{get_cross_section(nuclide_ID, reaction = 103)}\n")

    return

# Define the Element Symbols. Eisteinium included to prevent list index issues when calculating the parent nuclide name of Cf.
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']

# Reads in the csv file containing the data to be saved to the MPR file and assigns it to an array
#file_dir = r'C:\Users\sam.taylor\OneDrive - Newcleo\Documents\Modelling_LFR\Generating_MPR_file\ZAID_results.csv'
file_dir = r'/Users/sam/Documents/NewcleoInternship/ZAID_results.csv'
df = pd.read_csv(file_dir, header = 0)

# Assigns each column of csv file to a variable
nuclide = df[df.columns[0]]
zaid = df[df.columns[1]]
orion_id = df[df.columns[2]]
NumberReactions = df[df.columns[3]]
daughter_16 = df[df.columns[4]]
daughter_1