# This file should loop through all the endf data and create a csv file with:
# column 1 = endf ID
# column 2 = 



# For parent IDs, calculate the possible parent IDs for all reactions of a given nuclide then check if these IDs are in the database

import numpy as np
import endf
import os

# Define the Element Symbols
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf']

def build_nuclide_ID(atomic_number, mass_number):
    """_summary_

    Args:
        atomic_number (_type_): _description_
        mass_number (_type_): _description_

    Returns:
        _type_: _description_
    """
    if A < 10:
        ID = str(atomic_number) + "00" + str(mass_number)
    else:
        if A < 100:
            ID = str(atomic_number) + "0" + str(mass_number)
        else:
            ID = str(atomic_number) + str(mass_number)
    
    # Nuclear IDs in form ZZAAAI, where I = 0,1,2. I > 0 is for metastable istotopes. For now ignoring these as not in ENDF database.
    ID = str(ID) + "0"
    
    return ID

def calculate_daughter(atomic_number, mass_number, reaction):
    """_summary_

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
        return daughter_ID_16
    
    elif reaction ==17:
        # If nuclide is parent of reaction 17 (n,3n) then daughter should have A-2
        daughter_A_17 = mass_number - 2
        daughter_ID_17 = build_nuclide_ID(atomic_number, daughter_A_17)
        return daughter_ID_17
    
    elif reaction == 18:
        # If nuclide is parent of reaction 18 (n,fission) then daughter ID is 0
        daughter_ID_18 = 0
        return daughter_ID_18
    
    elif reaction == 102:
        # If nuclide is parent of reaction 102 then daughter should have A+1
        daughter_A_102 = mass_number + 1
        daughter_ID_102 = build_nuclide_ID(atomic_number, daughter_A_102)
        return daughter_ID_102
    
    elif reaction == 103:
        # If nuclide is parent of reaction 103 then daughter should have same A and Z-1
        daughter_Z_103 = atomic_number - 1
        daughter_ID_103 = build_nuclide_ID(daughter_Z_103, mass_number)
        return daughter_ID_103

    return 


# Path to Lib80x files, which are the cross-sections of reactions in the ENDF database
lib_dir = r'C:\USers\sam.taylor\OneDrive - Newcleo\Documents\Modelling_LFR\Generating_MPR_file\Lib80x\Lib80x'

nuclide = 1001
while nuclide < 98254:
    nuclide_str = str(nuclide)

    # Set Z to the first digit if the second digit is '0' and the nuclide ID number is less than 9999
    if nuclide < 10000:
        Z = int(nuclide_str[0])
    else:
        # Otherwise, set Z to the first two digits
        Z = int(nuclide_str[:2])
    
    # Mass number is set to figure left after removing first 2 digits
    A = nuclide_str[2:]

    # Remove leading zeros from A; if A is empty after stripping, set it to 0
    A = int(A.lstrip('0')) if A.lstrip('0') else 0
    
    # Nuclide ID in form ZZAAAI
    nuclide_ID = build_nuclide_ID(Z, A)
    
    # Build infile name based on chosen nuclide
    infile_str = (lib_dir + element_symbols[Z - 1] + "/" + nuclide_str + ".802nc")
    
    # Check whether the nuclide is in the database at all by validating if ENDF file exists
    if os.path.isfile(infile_str):
        table = endf.ace.get_table(infile_str)

        # Check whether the database has the required reactions for the nuclide and if so calculate the daughter nuclide ID in the form ZZAAAI
        if 16 in (table.interpret().reactions).keys():
            daughter_16 = calculate_daughter(Z, A, reaction = 16)
        if 17 in (table.interpret().reactions).keys():
            daughter_17 = calculate_daughter(Z, A, reaction = 17)
        if 18 in (table.interpret().reactions).keys():
            daughter_18 = calculate_daughter(Z, A, raction = 18)
        if 102 in (table.interpret().reactions).keys():
            daughter_102 = calculate_daughter(Z, A, reaction = 102)
        if 103 in (table.interpret().reactions).keys():
            daughter_103 = calculate_daughter(Z, A, reaction = 103)
    

        # If nuclide is daughter of reaction 16 (n,2n) then parent should have A+1
        parent_A_16 = A + 1
        parent_16 = build_nuclide_ID(Z, parent_A_16)
        
        # If nuclide is daughter of reaction 17 (n,3n) then parent should have A+2
        parent_A_17 = A + 2
        parent_17 = build_nuclide_ID(Z, parent_A_17)
        
        # If nuclide is daughter of reaction 18 (n,fission) then parent ID is 0
        parent_18 = 0
        
        # If nuclide is daughter of reaction 102 then parent should have A-1
        parent_A_102 = A - 1
        parent_102 = build_nuclide_ID(Z, parent_A_102)
        
        # If nuclide is daughter of reaction 103 then parent should have same A and Z+1
        parent_Z_103 = Z + 1
        parent_103 = build_nuclide_ID(parent_Z_103, A)
        
        
        
        
        
        
        
    nuclide = nuclide + 1