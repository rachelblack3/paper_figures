# Code that compares the power in survey and burst 

# make 12x12 plot

# numerical essentials 
import numpy as np

# for plotting
from matplotlib import pyplot as plt, animation
from metpy.plots import add_timestamp
import os
import pandas as pd
import xarray as xr

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
# for cdf reading
from spacepy import toolbox
#spacepy.toolbox.update(leapsecs=True)
from spacepy import pycdf
import h5py

# other numerical tools
import os
import glob
from datetime import datetime,timedelta,date
import string

plt.style.use("ggplot")

def to_dt(nanoseconds_value):
    # the default time format in datasets
    # Convert nanoseconds to seconds
    seconds_value = nanoseconds_value / 1e9

    # Unix epoch
    epoch = datetime(1970, 1, 1)

    # Convert to datetime
    human_readable_datetime = epoch + timedelta(seconds=seconds_value)

    return human_readable_datetime

def make_stat_histogram(stat,variable1,variable2):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    N=10
    AE_bins = np.nanquantile(variable1, np.linspace(0, 1, N + 1))

    fpe_fce_bins = np.nanquantile(variable2, np.linspace(0, 1, N + 1))

    power = stat

    # Digitize the MLT and L values to get the bin indices
    AE_bin_indices = np.digitize(variable1, AE_bins) -1 # Subtract 1 to make bins zero-indexed

    fpe_fce_bin_indices = np.digitize(variable2, fpe_fce_bins) -1
    

    # Create a 2D array to store the binned power data
    binned_power = np.zeros((len(AE_bins)-1, len(fpe_fce_bins)-1))
    bin_counts = np.zeros((len(AE_bins)-1, len(fpe_fce_bins)-1))

    # Aggregate power values by bin
    for i in range(len(power)):
        ae_idx = AE_bin_indices[i]
        fpe_fce_idx = fpe_fce_bin_indices[i]
        if 0 <= ae_idx < len(AE_bins)-1 and 0 <= fpe_fce_idx < len(fpe_fce_bins)-1:  # Ensure valid indices
            if np.isnan(power[i])==False:
                binned_power[ae_idx, fpe_fce_idx] += power[i]
                bin_counts[ae_idx, fpe_fce_idx] += 1

    # Optionally, calculate the mean by counting occurrences in each bin
    #bin_counts = np.histogram2d(variable1, variable2, bins=[AE_bins, fpe_fce_bins])[0]
    #for count in bin_counts:
        #print(count)
    binned_mean_power = np.divide(binned_power, bin_counts, where=bin_counts > 50)  # Avoid division by zero
    # To replace the result with NaN where the condition isn't met:
    binned_mean_power = np.where(bin_counts > 50, binned_mean_power, np.nan)

    # 'binned_mean_power' now holds the average power in each MLT-L bin


    return binned_mean_power,AE_bins,fpe_fce_bins

# making a plot of:

# - each band (3 - low frequency, lowerband, upperband)
# - and each stat (2 - IQR_norm and ratio)
# - and each latitude (2 - equatorial and off-equatorial)

# Lets have 3 columns, 2 rows

# Concatenate along the 'x' dimension
chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24
# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"]
alphabet = list(string.ascii_lowercase)

MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)
# create condition figure
fig = plt.figure(figsize=(21, 14))
# GridSpec: 3 rows, 2 columns
gs = GridSpec(1, 2, figure=fig, width_ratios=[0.5,8])


name = 'lowerband'
# Top row: 3 subfigures, bottom is 2
subfig = fig.add_subfigure(gs[1])

#slice_fig = fig.add_subfigure(gs[3])
axs = subfig.subplots(1, 1,)

# set conditions 
upperband_cond = (np.isnan(chorus_result["upperband_IQR"]/chorus_result["upperband_median"])==False)
conditions = (chorus_result[f'lowerband_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_survey_power'] !=0)& (chorus_result[f'upperband_peak']<1e4)&(chorus_result[f'upperband_peak'] !=0)& MLAT_cond & upperband_cond

# Create a figure
# Create the figure
# Create a figure and axes array for subplots


iqr_vals_upper =  np.where(conditions,chorus_result[f'upperband_IQR'][:].values, np.nan)
median_vals_upper = np.where(conditions,chorus_result[f'lowerband_median'][:].values, np.nan)
iqr_norm_upper = iqr_vals_upper/median_vals_upper
iqr_norm_upper = iqr_norm_upper



iqr_vals_lower =  np.where(conditions,chorus_result[f'lowerband_IQR'][:].values, np.nan)
median_vals_lower = np.where(conditions,chorus_result[f'lowerband_median'][:].values, np.nan)
iqr_norm_lower = iqr_vals_lower/median_vals_lower
iqr_norm_lower = iqr_norm_lower

for i,val in enumerate(iqr_norm_upper):
    if val >1e4:
        iqr_norm_upper[i] = np.nan
        iqr_norm_lower[i] = np.nan
        
axs.scatter(iqr_norm_lower,iqr_norm_upper)




plt.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig7/upperversuslower.png')






