# Code that compares the power in survey and burst 

# make 12x12 plot

# numerical essentials 
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
plt.rcParams['image.cmap'] = 'plasma'

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
            if ~np.isnan(power[i]) and power[i]<1e4:
                binned_power[ae_idx, fpe_fce_idx] += power[i]
                bin_counts[ae_idx, fpe_fce_idx] += 1

    print("bin counts are (ae,fpefce)",bin_counts)

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
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')
burstA = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/A/full_12_18_withTotalPowerAndSurv.nc')
burstB = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/B/full_12_19_withTotalPowerAndSurv.nc')
chorus_result = xr.concat([burstA, burstB], dim="x")
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24
# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"]
alphabet = list(string.ascii_lowercase)

MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)
# create condition figure
fig = plt.figure(figsize=(22, 12))
# GridSpec: 3 rows, 2 columns
gs = GridSpec(1, 7, figure=fig, width_ratios=[0.8,2,0.3,2,0.3,2,0.8])
print(f'{names[2]}')

flags = ['lf','lb','ub']

for m, (i,name,flag) in enumerate(zip([1,3,5],names,flags)):

    # Top row: 3 subfigures, bottom is 2
    subfig = fig.add_subfigure(gs[i])

    #slice_fig = fig.add_subfigure(gs[3])
    axs = subfig.subplots(2, 1,)
    subfig.subplots_adjust(wspace=0.4, hspace=0.3)
    axsTop = axs[0]
    axsBottom = axs[1]
    #axsSlice = slice_fig.subplots(1, 1,) 

    axsTop.set_title(f'{alphabet[m]})', loc='left', fontsize=20)
    axsBottom.set_title(f'{alphabet[m+3]})', loc='left', fontsize=20)

    print(f'{name} started')
    # set conditions 
    conditions = (chorus_result[f'{flag}_flag']==1) & MLAT_cond & (chorus_result[f'{name}_peak']<1e2)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_surv_power'] !=0)

    # Create a figure
    # Create the figure
    # Create a figure and axes array for subplots


    iqr_vals =  np.where(conditions,chorus_result[f'{name}_iqr'][:].values, np.nan)
    median_vals = np.where(conditions,chorus_result[f'{name}_med'][:].values, np.nan)
    survey_vals = np.where(conditions,chorus_result[f'{name}_surv_power'][:].values, np.nan)
    peak_vals = np.where(conditions,chorus_result[f'{name}_peak'][:].values, np.nan)
    times = np.where(conditions,chorus_result[f'timestamp'][:].values, np.nan)
    iqr_over_median = iqr_vals/median_vals
    ratio_vals = peak_vals/survey_vals

    
    if name=="lowerband":
        for j in range(50):
            if ~np.isnan(peak_vals[j]):
                print(j,"sanity check:",peak_vals[j],survey_vals[j],ratio_vals[j],to_dt(times[j]))

    if name=="upperband":
        for j in range(50):
            if ~np.isnan(peak_vals[j]):
                print(j,"sanity check:",peak_vals[j],survey_vals[j],ratio_vals[j])


    # do the ratio on the top plots
    
    stat = ratio_vals
    fpe_fce_vals = np.where(conditions&(chorus_result[f'fpe_fce']<20),chorus_result[f'fpe_fce'][:].values, np.nan)
    ae_vals = np.where(conditions,chorus_result[f'AE'][:].values, np.nan)

    # do the left figure
    title = r'$P_{burst,max}/P_{survey,FFT}$'
    ae_fpe_array, ae_bins, fpe_fce_bins = make_stat_histogram(stat,ae_vals,fpe_fce_vals)
    ratio_norm=mcolors.LogNorm(vmin=1e0,vmax=5*1e2)
    c = axsTop.pcolor(ae_fpe_array,norm=ratio_norm)

    if i ==1:
        axsTop.set_ylabel("AE",fontsize=20)

        # Add side text to Plot 1
        axsTop.text(-0.42, 0.5,                      # x, y in axes coordinates
        r'$P_{burst,max}/P_{survey,FFT}$',
        transform=axsTop.transAxes,
        rotation=90,
        va='center',
        ha='left',
        fontsize=20,
        color='black'
        )

        axsTop.text(0., 1.2,                      # x, y in axes coordinates
        r'low frequency',
        transform=axsTop.transAxes,
        rotation=0,
        va='center',
        ha='left',
        fontsize=20,
        color='black'
        )

        axsBottom.text(-0.42, 0.5,                      # x, y in axes coordinates
        r'$IQR_{norm}$',
        transform=axsBottom.transAxes,
        rotation=90,
        va='center',
        ha='left',
        fontsize=20,
        color='black'
        )

    else:
        axsTop.text(0., 1.2,                      # x, y in axes coordinates
        f'{name}',
        transform=axsTop.transAxes,
        rotation=0,
        va='center',
        ha='left',
        fontsize=20,
        color='black'
        )

        
    #axsTop.set_title(f'{alphabet[i]} i.', loc='left', fontsize=14)
    #axsBottom.set_title(f'ii.', loc='left', fontsize=14)


    # Calculate bin centers
    ae_centers = 0.5 * (ae_bins[:-1] + ae_bins[1:])
    fpe_fce_centers = 0.5 * (fpe_fce_bins[:-1] + fpe_fce_bins[1:])

    # Set ticks at the center of each bin
    axsTop.set_yticks(np.arange(len(ae_bins)))
    axsTop.set_yticklabels([f"{val:.1f}" for val in ae_bins], fontsize=20)

    axsTop.set_xticks(np.arange(len(fpe_fce_bins)))
    axsTop.set_xticklabels([f"{val:.2f}" for val in fpe_fce_bins],rotation=45, fontsize=20)

    # do the IQR on the bottom plots
    stat = iqr_over_median
    fpe_fce_vals = np.where(conditions&(~np.isnan(iqr_over_median)&(chorus_result[f'fpe_fce']<20)),chorus_result[f'fpe_fce'][:].values, np.nan)
    ae_vals = np.where(conditions&(~np.isnan(iqr_over_median)),chorus_result[f'AE'][:].values, np.nan)

    cbar_labels = [r'$nT^2$',r'$nT^2$',r'$Ratio$',r'$Count$']
    # do the left figure
    title = r'$IQR_{norm}$'
    ae_fpe_array, ae_bins, fpe_fce_bins = make_stat_histogram(stat,ae_vals,fpe_fce_vals)
    norm=mcolors.LogNorm(vmin=5*1e-1,vmax=1e2)
    c = axsBottom.pcolor(ae_fpe_array,norm=norm)

    if i ==1 or i==3:
        axsTop.set_ylabel("AE",fontsize=20)
        axsBottom.set_ylabel("AE",fontsize=20)
    
    axsBottom.set_xlabel(r'$f_{pe}/f_{ce}$',fontsize=20)

    # Calculate bin centers
    ae_centers = 0.5 * (ae_bins[:-1] + ae_bins[1:])
    fpe_fce_centers = 0.5 * (fpe_fce_bins[:-1] + fpe_fce_bins[1:])

    # Set ticks at the center of each bin
    axsBottom.set_yticks(np.arange(len(ae_bins)))
    axsBottom.set_yticklabels([f"{val:.1f}" for val in ae_bins],fontsize=20)
    axsBottom.set_xticks(np.arange(len(fpe_fce_bins)))
    axsBottom.set_xticklabels([f"{val:.2f}" for val in fpe_fce_bins],rotation=45, fontsize=20)

    if i ==5:
        print("i am here")
        # Add a separate axis for the color bar of the IQR plot
        sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
        sm.set_array([])  # Required for the colorbar to work properly
        cbar_ax = subfig.add_axes([0.95, 0.175, 0.05, 0.2])  # [left, bottom, width, height]
        # Make the axes invisible
        for spine in cbar_ax.spines.values():
            spine.set_visible(False)
        cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        # Add the color bar to the figure
        #cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
        cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
        cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
        cbar.set_label(r'$IQR_{norm}$', fontsize=20) 
        cbar_ax.set_facecolor('none') 

        # Add a separate axis for the color bar of the IQR plot
        sm = plt.cm.ScalarMappable(cmap='plasma', norm=ratio_norm)
        sm.set_array([])  # Required for the colorbar to work properly
        cbar_ax = subfig.add_axes([0.95, 0.6, 0.05, 0.2])  # [left, bottom, width, height]
        # Make the axes invisible
        for spine in cbar_ax.spines.values():
            spine.set_visible(False)
        cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        # Add the color bar to the figure
        #cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
        cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
        cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
        cbar.set_label(r'$P_{burst,max}/P_{survey,FFT}$', fontsize=20) 
        cbar_ax.set_facecolor('none') 


    plt.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig7/NewStats/low_6x6BIGGERV1_1.png')
