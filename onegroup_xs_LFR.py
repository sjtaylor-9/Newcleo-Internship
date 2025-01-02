"""_summary_
"""
# ---------  Import Libraries  --------- #
import numpy as np
from numpy import genfromtxt
import endf
import pandas as pd
import os
import argparse
# --------- Global Variables ----------- #
element_symbols = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',	'Pd','Ag','Cd','In','Sn','Sb','Te',	'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No']
# Set up units (working in CGS system)
barn = 1.e-24 #cm^2
# MT IDs of reactions under consideration
reaction_ids = [16, 17, 18, 102, 103]
reaction_strs = ["(n,2n)", "(n,3n)", "(n,fission)", "(n,$\gamma$)", "(n,p)"]
# -------------  Functions  ------------ #
def parse_arguments():
    """
    Parses the arguments needed along the code. Arguments:

    --reactor               Used to specify which reactor out of the LFR200 and LFR30 is being examined.
                            The argument must be one of: [LFR30, LFR200].
    --enrichment_zone       If the LFR200 is being examined, then the enrichment zone must be specified. The data for the LFR30 neutron spectrum only describes one enrichment zone, the fuel assembly.
                            The argument must be one of: [inner, middle, outer, FA].
    --path                  Used to specify the directory in which the one-group cross-section data should be written. It is not required,
                            in the case it is not specified, the default path is the current working directory.
    --newcleo_input         Used to specify the directory in which internal newcleo inpput data should be found. This includes data such as the neutron spectra.
    --cross_section_data    Used to specify the directory in which the public data regarding the ENDF and JEFF reaction cross-sections can be found.
    --reaction_data
    --burnup            
    
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
        "--newcleo_input",
        type=dir_path,
        required=True,
        help="flag to set the path where internal newcleo input data should be found"
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

def get_neutron_spectrum():
    
    if args.reactor == 'LFR30':
        # csv file with the edges of the energy groups
        neutronspec_bin_edges = genfromtxt(f'{args.newcleo_input}/LFR30_neutronspec_groupedges.csv', delimiter = ',')
        # group cross-section data according to the neutron spec group structure
        central_value_neutron_flux = genfromtxt(f'{args.newcleo_input}/LFR30_neutronspec_values.csv', delimiter = ',')
    elif args.reactor == 'LFR200':
        if args.enrichment_zone == 'inner':
            neutronspec_bin_edges = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'INNER FUEL',
                          usecols = ['Emax'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy() # Converts pd dataframe to numpy array
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1) # Ensures the numpy array is 1D
            central_value_neutron_flux = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'INNER FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()
        elif args.enrichment_zone == 'middle':
            neutronspec_bin_edges = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'MIDDLE FUEL',
                          usecols = ['Emax'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy()
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1)
            central_value_neutron_flux = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'MIDDLE FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()
        elif args.enrichment_zone == 'outer':
            neutronspec_bin_edges = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = ['Emax'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy()
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1)
            central_value_neutron_flux = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()

    num_energy_bins = len(neutronspec_bin_edges) - 1
    return neutronspec_bin_edges, num_energy_bins, central_value_neutron_flux

def get_jeff_reaction_xs(bin_edges, num_bins, neutron_flux):
    """
    Get JEFF 3.3 reaction cross-sections
    """
    # Select element, isotope and reaction (by MT ID)
    lib_dir = f'{args.cross_section_data}/ace_900/ace_900/'
    jeff_reactions = genfromtxt(f'{args.reaction_data}/JEFF_reactions.csv', delimiter = ',')
    n_reactions = len(jeff_reactions[:, 0])
    jeff_cross_sections = []

    # Get Z from element symbol
    for i in range(n_reactions):

        parent_ZAI = jeff_reactions[i, 0]

        test_reaction_id = jeff_reactions[i, 1]

        meta_flag = np.mod(parent_ZAI, 10)
        #print(meta_flag)
        #print("%i, %s, %i" % (parent_ZAI,ZAI_str,meta_flag))
        parent_ZA = parent_ZAI - meta_flag

        parent_Z = int(parent_ZA / 10000)
        parent_A = int(((parent_ZA) - (10000 * parent_Z)) / 10)

        if(meta_flag == 1): meta_str = "m"
        else: meta_str = "g"

        element = element_symbols[parent_Z - 1]
        filename = "%i-%s-%i%s-900.ace" % (parent_Z, element, parent_A, meta_str)
        infile_str = lib_dir + filename
        table = endf.ace.get_table(infile_str)
        test_reaction = (table.interpret().reactions[test_reaction_id].xs)['900K']

        d1 = {'Energy': test_reaction.x, 'Sigma': test_reaction.y}
        test_df = pd.DataFrame(data = d1)

        # group cross-section data according to the neutron spec group structure
        grouped_xs_df = test_df.groupby(pd.cut(test_df["Energy"], bin_edges), observed = False).mean()
        grouped_xs_df = grouped_xs_df.Sigma.reset_index()
        grouped_xs = grouped_xs_df["Sigma"].copy()

        # bit of a bodge - fixes any empty bins in the middle of the range
        for i in range(1, num_bins - 1):
            if np.isnan(grouped_xs[i]):
                grouped_xs[i] = 0.5 * (grouped_xs[i - 1] + grouped_xs[i + 1])

        prod_sum = 0
        flux_sum = 0
        for i in range(len(grouped_xs)):
            if(np.isnan(grouped_xs[i])): grouped_xs[i] = 0
            prod_sum += grouped_xs[i] * neutron_flux[i]#*neutronspec_binwidths[i]
            flux_sum += neutron_flux[i]#*neutronspec_binwidths[i]


        if test_reaction_id == 16:
            daughter_Z = parent_Z
            daughter_A = parent_A - 1
        if test_reaction_id == 17:
            daughter_Z = parent_Z
            daughter_A = parent_A - 2
        if test_reaction_id == 102:
            daughter_Z = parent_Z
            daughter_A = parent_A + 1
        if test_reaction_id == 103:
            daughter_Z = parent_Z - 1
            daughter_A = parent_A
        if test_reaction_id == 18:
            daughter_Z = 0
            daughter_A = 0
        daughter_ZAI = (daughter_Z * 10000) + (daughter_A * 10)
        #print("%i,%i,%i,%.5e" % (parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)))
        jeff_cross_sections.append([parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)])

    return jeff_cross_sections

def get_endf_reaction_xs(bin_edges, num_bins, neutron_flux):
    """
    Get remaining ENDF reaction cross-sections
    """
    # Select element, isotope and reaction (by MT ID)
    lib_dir = f'{args.cross_section_data}/Lib80x/Lib80x/'
    missing_reactions = genfromtxt(f'{args.reaction_data}/reactions_in_endf_not_jeff.csv', delimiter = ',')
    n_reactions = len(missing_reactions[:, 0])
    endf_cross_sections = []
    
    # Get Z from element symbol
    for i in range(n_reactions):
        parent_ZAI = missing_reactions[i, 0]

        test_reaction_id = missing_reactions[i, 1]

        meta_flag = np.mod(parent_ZAI, 10)
        #print("%i, %s, %i" % (parent_ZAI,ZAI_str,meta_flag))
        file_ZA = int(parent_ZAI/ 10)
        parent_Z = int(parent_ZAI/ 10000)
        parent_A = int((parent_ZAI - (10000 * parent_Z))/ 10)

        if(meta_flag == 1): meta_str = "m1_"
        else: meta_str = ""

        element = element_symbols[parent_Z - 1]
        #print(meta_str)
        #print(file_ZA)

        filename = "%s/%s%i.802nc" % (element, meta_str, file_ZA)
        infile_str = lib_dir + filename
        table = endf.ace.get_table(infile_str)
        test_reaction = (table.interpret().reactions[test_reaction_id].xs)['900K']

        d1 = {'Energy': test_reaction.x, 'Sigma': test_reaction.y}
        test_df = pd.DataFrame(data = d1)

        # group cross-section data according to the neutron spec group structure
        grouped_xs_df = test_df.groupby(pd.cut(test_df["Energy"], bin_edges), observed = False).mean()

        grouped_xs_df = grouped_xs_df.Sigma.reset_index()
        grouped_xs = grouped_xs_df["Sigma"].copy()

        # bit of a bodge - fixes any empty bins in the middle of the range
        for i in range(1, num_bins - 1):
            if np.isnan(grouped_xs[i]):
                grouped_xs[i] = 0.5 * (grouped_xs[i - 1] + grouped_xs[i + 1])


        prod_sum = 0
        flux_sum = 0
        for i in range(len(grouped_xs)):
            if(np.isnan(grouped_xs[i])): grouped_xs[i] = 0
            prod_sum += grouped_xs[i] * neutron_flux[i]#*neutronspec_binwidths[i]
            flux_sum += neutron_flux[i]#*neutronspec_binwidths[i]


        if test_reaction_id == 16:
            daughter_Z = parent_Z
            daughter_A = parent_A - 1
        if test_reaction_id == 17:
            daughter_Z = parent_Z
            daughter_A = parent_A - 2
        if test_reaction_id == 102:
            daughter_Z = parent_Z
            daughter_A = parent_A + 1
        if test_reaction_id == 103:
            daughter_Z = parent_Z - 1
            daughter_A = parent_A
        if test_reaction_id == 18:
            daughter_Z = 0
            daughter_A = 0
        daughter_ZAI = (daughter_Z * 10000) + (daughter_A * 10)
        #print("For reaction %i on nuclide %s-%i, onegroup XS = %.3e" % (test_reaction_id,element,test_A,(prod_sum/flux_sum)))
        #print("%s,%i,%i,%i,%s,%i,%.5e" % (element,test_Z,test_A,metastable_flag,reaction_strs[reaction_ids.index(test_reaction_id)],test_reaction_id,(prod_sum/flux_sum)))
        #print("%i,%i,%i,%.5e" % (parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)))
        endf_cross_sections.append([parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)])

    return endf_cross_sections

def save_results(jeff_xs, endf_xs):
    """_summary_

    Args:
        jeff_xs (_type_): _description_
        endf_xs (_type_): _description_
    """
    # Convert the cross-sections to standard form with 2 d.p.
    for i in range(len(endf_xs)):
        xs = endf_xs[i][3]
        formatted_xs = np.format_float_scientific(xs, precision = 2)
        endf_xs[i][3] = formatted_xs
    for i in range(len(jeff_xs)):
        xs = jeff_xs[i][3]
        formatted_xs = np.format_float_scientific(xs, precision = 2)
        jeff_xs[i][3] = formatted_xs
    # Convert numpy arrays to pandas dataframes
    jeff_data_df = pd.DataFrame(jeff_xs)
    endf_data_df = pd.DataFrame(endf_xs)

    # Format the dataframes such that first 2 columns don't have .0 decimals
    jeff_data_df[0] = jeff_data_df[0].apply(lambda x: f'{x: .0f}')
    jeff_data_df[1] = jeff_data_df[1].apply(lambda x: f'{x: .0f}')
    endf_data_df[0] = endf_data_df[0].apply(lambda x: f'{x: .0f}')
    endf_data_df[1] = endf_data_df[1].apply(lambda x: f'{x: .0f}')
    
    # Create directories for output files
    jeff_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_JEFF_reactiondata.csv'
    endf_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_ENDF_reactiondata.csv'

    # Outputs the JEFF and ENDF reaction cross-sections to .csv files
    jeff_data_df.to_csv(jeff_file, index = False, header = False)
    endf_data_df.to_csv(endf_file, index = False, header = False)

    # Creates a single numpy array with both the JEFF and ENDF reaction cross-sections and converts this to a pandas dataframe
    all_reactions = np.vstack((jeff_xs, endf_xs))
    all_reactions_df = pd.DataFrame(all_reactions)

    # Sort the dataframe so that it is in the order of ascending ZAIDs
    all_reactions_df = all_reactions_df.sort_values(by = all_reactions_df.columns[0])
    all_reactions_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_all_reactiondata.csv'
    all_reactions_df.to_csv(all_reactions_file, index = False, header = False)

    return
# -------------  Main Code  ------------ #
args = parse_arguments()

neutronspec_edges, n_energy_bins, neutronspec_central = get_neutron_spectrum()
print(f"Calculating the cross-sectional data for the JEFF reactions")
jeff_data = get_jeff_reaction_xs(neutronspec_edges, n_energy_bins, neutronspec_central)
print(f"Calculating the cross-sectional data for the ENDF reactions")
endf_data = get_endf_reaction_xs(neutronspec_edges, n_energy_bins, neutronspec_central)

save_results(jeff_data, endf_data)
print(f"Calculated the cross-sectional data for the {args.reactor} in the enrichment zone: {args.enrichment_zone} at the {args.burnup} burnup step")