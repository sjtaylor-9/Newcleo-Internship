# ---------  Import Libraries  --------- #
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.serif'] = ['Times New Roman']
import numpy as np
# -------------  Functions  ------------ #
def cdf_total_fuel_mass(dataframe, first_reactor_deployed):
    """_summary_

    Args:
        dataframe (_type_): _description_
        first_reactor_deployed (_type_): _description_
    """
    # These parameters are allowed to vary to test supplement business plan
    uk_Pu_stockpile = 140
    last_deployment_year = 2050
    reactor_lifetime = 60
    first_shutdown = first_reactor_deployed + reactor_lifetime
    
    # Calculate the cumulative total Pu mass
    dataframe['Cumulative Total Pu'] = dataframe['Total Pu'].cumsum()
    
    # Plots the cumulative total as function of data
    fig = plt.figure(figsize = (10, 6))
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu'],
             color = 'black')
    
    # Interpolates the year for which the UK Pu stockpile is used up in
    interpolated_year_stockpile_used = np.interp(uk_Pu_stockpile, df['Cumulative Total Pu'], df['Date'])
    print('The UK Plutonium stockpile will be exhausted in', interpolated_year_stockpile_used)
    # Plot the interpolated point and add lines to the axes to show when this occurs
    plt.scatter(interpolated_year_stockpile_used,
                uk_Pu_stockpile,
                color = 'skyblue',
                zorder = 5,
                label = f'UK Pu stockpile ({uk_Pu_stockpile} t)')
    # plt.hlines(y = uk_Pu_stockpile,
    #             xmin = dataframe['Date'].min(),
    #             xmax = interpolated_year_stockpile_used,
    #             colors = 'orange',
    #             linestyles = '--',
    #             linewidth = 1.4)
    # plt.vlines(x = interpolated_year_stockpile_used,
    #             ymin = dataframe['Cumulative Total Pu'].min(),
    #             ymax = uk_Pu_stockpile,
    #             colors = 'orange',
    #             linestyles = '--',
    #             linewidth = 1.4)
    
    # Plot the point at which all reactors have been deployed
    plt.scatter(last_deployment_year,
                dataframe[dataframe['Date'] == last_deployment_year]['Cumulative Total Pu'].values[0],
                color = 'red',
                zorder = 5,
                label = 'Full fleet deployed')
    
    # Plot the point at which the first set of reactors have been shut down
    plt.scatter(first_shutdown,
                dataframe[dataframe['Date'] == first_shutdown]['Cumulative Total Pu'].values[0],
                color = 'limegreen',
                zorder = 5,
                label = 'Beginning of fleet shutdown')
    
    # Define axes limits
    plt.xlim(first_reactor_deployed, dataframe['Date'].max())
    plt.ylim(0, dataframe['Cumulative Total Pu'].max() + 100)
    
    # Axes and legend labels
    plt.title('Cumulative total Pu mass used in LFR200 fleet')
    plt.xlabel('Date [yr]')
    plt.ylabel('Cumulative total Pu mass [t]')
    plt.legend(loc = 'upper left')
    
    file_path = r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/LFR200_Simulation/Testing_Reactor_Fleet/Images/cumulative_Pu_mass_fleet.png'
    plt.savefig(file_path, bbox_inches = "tight", dpi = 300)
    return
def cdf_individual_fuel_masses(dataframe, first_reactor_deployed):
    
    # Calculate the cumulative total mass of each Pu isotope and Am241
    dataframe['Cumulative Total Pu238'] = dataframe['3+PU238'].cumsum()
    dataframe['Cumulative Total Pu239'] = dataframe['3+PU239'].cumsum()
    dataframe['Cumulative Total Pu240'] = dataframe['3+PU240'].cumsum()
    dataframe['Cumulative Total Pu241'] = dataframe['3+PU241'].cumsum()
    dataframe['Cumulative Total Pu242'] = dataframe['3+PU242'].cumsum()
    dataframe['Cumulative Total Am241'] = dataframe['3+AM241'].cumsum()

    # Plots the cumulative total as function of data
    fig = plt.figure(figsize = (10, 6))
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu238'],
             color = 'black',
             label = 'Pu238')
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu239'],
             color = 'skyblue',
             label = 'Pu239')
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu240'],
             color = 'red',
             label = 'Pu240')
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu241'],
             color = 'green',
             label = 'Pu241')
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu242'],
             color = 'purple',
             label = 'Pu242')
    plt.plot(dataframe['Date'],
             dataframe['Cumulative Total Pu238'],
             color = 'orange',
             label = 'Am241')

    # Define axes limits
    plt.xlim(first_reactor_deployed, dataframe['Date'].max())
    
    # Axes and legend labels
    plt.title('Cumulative total fuel mass of Pu vector in LFR200 fleet')
    plt.xlabel('Date [yr]')
    plt.ylabel('Cumulative total mass of each nuclide in Pu vector [t]')
    plt.legend(loc = 'upper left')
    
    file_path = r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/LFR200_Simulation/Testing_Reactor_Fleet/Images/cumulative_Pu_mass_ind_nuclides_fleet.png'
    plt.savefig(file_path, bbox_inches = "tight", dpi = 300)
    
    # Create a second plot with a logy axis
    plt.yscale('log')
    plt.ylabel('log(Cumulative total mass of each nuclide in Pu vector [t])')
    file_path = r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/LFR200_Simulation/Testing_Reactor_Fleet/Images/cumulative_Pu_mass_ind_nuclides_fleet_logplot.png'
    plt.savefig(file_path, bbox_inches = "tight", dpi = 300)
    return
# -------------  Main Code  ------------ #
# Read in the .xlsx file containing the refuelling dates for each reactor and assign it to a pandas dataframe
file_dir = r'/mnt/c/Users/sam.taylor/OneDrive - Newcleo/Documents/Modelling_LFR/LFR200_Simulation/Testing_Reactor_Fleet/Pu_fuel_mass_per_year.csv'
df = pd.read_csv(file_dir, header = 0)

# Set these as variables so can easil change if deployment schedule changes
first_deployment_year = 2033
fuel_fabrication_time = 2

# Convert the refuel dates into relevant years
df['Date'] = df['Date'] + first_deployment_year - fuel_fabrication_time
# Calculate the total mass of Pu at each timestep
df['Total Pu'] = df['3+PU238'] + df['3+PU239'] + df['3+PU240'] + df['3+PU241'] + df['3+PU242']

cdf_total_fuel_mass(df, first_deployment_year)
cdf_individual_fuel_masses(df, first_deployment_year)