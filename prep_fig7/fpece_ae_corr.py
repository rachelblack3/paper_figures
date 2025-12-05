# Code that compares the power in survey and burst 

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

plt.style.use("ggplot")

def make_stat_histogram(stat,variable1,variable2):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    N=10
    AE_bins = np.nanquantile(variable1, np.linspace(0, 1, N + 1))

    fpe_fce_bins = np.nanquantile(variable2, np.linspace(0, 1, N + 1))
    print(fpe_fce_bins,AE_bins)

    power = stat

    # Digitize the MLT and L values to get the bin indices
    AE_bin_indices = np.digitize(variable1, AE_bins) -1 # Subtract 1 to make bins zero-indexed

    fpe_fce_bin_indices = np.digitize(variable2, fpe_fce_bins) -1
    

    # Create a 2D array to store the binned power data
    binned_power = np.zeros((len(AE_bins)-1, len(fpe_fce_bins)-1))

    # Aggregate power values by bin
    for i in range(len(power)):
        ae_idx = AE_bin_indices[i]
        fpe_fce_idx = fpe_fce_bin_indices[i]
        if 0 <= ae_idx < len(AE_bins)-1 and 0 <= fpe_fce_idx < len(fpe_fce_bins)-1:  # Ensure valid indices
            binned_power[ae_idx, fpe_fce_idx] += power[i]

    # Optionally, calculate the mean by counting occurrences in each bin
    bin_counts = np.histogram2d(variable1, variable2, bins=[AE_bins, fpe_fce_bins])[0]
    #for count in bin_counts:
        #print(count)
    binned_mean_power = np.divide(binned_power, bin_counts, where=bin_counts > 50)  # Avoid division by zero
    # To replace the result with NaN where the condition isn't met:
    binned_mean_power = np.where(bin_counts > 50, binned_mean_power, np.nan)

    # 'binned_mean_power' now holds the average power in each MLT-L bin


    return binned_mean_power,AE_bins,fpe_fce_bins

# Concatenate along the 'x' dimension
chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24
# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"]

MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)

print(f'{names[2]}')
for name in names:

    # create condition figure
    fig = plt.figure(figsize=(28, 12))
    # GridSpec: 1 rows, 4 columns
    gs = GridSpec(1, 3, figure=fig, width_ratios=[4,0.5,5])

    # Top row: 4 subfigures
    left_fig = fig.add_subfigure(gs[0])
    right_fig = fig.add_subfigure(gs[2])
    #slice_fig = fig.add_subfigure(gs[3])
    axsLeft = left_fig.subplots(1, 1,)
    axsRight = right_fig.subplots(2, 2,) 
    #axsSlice = slice_fig.subplots(1, 1,) 

    print(f'{name} started')
    # set conditions 
    conditions = (chorus_result[f'{name}_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_survey_power'] !=0) & MLAT_cond

    # Create a figure
    # Create the figure
    # Create a figure and axes array for subplots


    iqr_vals =  np.where(conditions,chorus_result[f'lowerband_IQR'][:].values, np.nan)
    median_vals = np.where(conditions,chorus_result[f'lowerband_median'][:].values, np.nan)
    ae_vals = np.where(conditions,chorus_result[f'AE'][:].values, np.nan)
    fpe_fce_vals = np.where(conditions,chorus_result[f'fpe_fce'][:].values, np.nan)
    iqr_over_median = iqr_vals/median_vals
    

    stat = iqr_over_median
    cbar_labels = [r'$nT^2$',r'$nT^2$',r'$Ratio$',r'$Count$']


    # do the left figure
    title = r'$IQR_{norm}$'
    ae_fpe_array, ae_bins, fpe_fce_bins = make_stat_histogram(stat,ae_vals,fpe_fce_vals)
    norm=mcolors.Normalize(vmin=0,vmax=30)
    c = axsLeft.pcolor(ae_fpe_array,norm=norm)
    axsLeft.set_title(title)
    axsLeft.set_ylabel("AE")
    axsLeft.set_xlabel(r'$f_{pe}/f_{ce}$')

    # Calculate bin centers
    ae_centers = 0.5 * (ae_bins[:-1] + ae_bins[1:])
    fpe_fce_centers = 0.5 * (fpe_fce_bins[:-1] + fpe_fce_bins[1:])

    # Set ticks at the center of each bin
    print(ae_bins)
    axsLeft.set_yticks(np.arange(len(ae_bins)))
    axsLeft.set_yticklabels([f"{val:.1f}" for val in ae_bins])

    axsLeft.set_xticks(np.arange(len(fpe_fce_bins)))
    axsLeft.set_xticklabels([f"{val:.2f}" for val in fpe_fce_bins])




    # Add a separate axis for the color bar of the IQR plot
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = left_fig.add_axes([0.95, 0.4, 0.05, 0.2])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$nT^2$', fontsize=14) 
    cbar_ax.set_facecolor('none') 

    # do the right figure

    # Plot histograms

    axsRight[1,0].hist(fpe_fce_vals[fpe_fce_vals<=20.],bins=10)

    axsRight[1,0].set_xlabel(r'$f_{pe}/f_{ce}$')

    axsRight[0,0].hist(ae_vals[ae_vals<=1500.],bins=10)

    axsRight[0,0].set_xlabel(r'$Ae\ index$')

    # do the quantile bins: 
    N = 10
    hist_counts, _ = np.histogram(fpe_fce_vals[fpe_fce_vals<=20.],bins=fpe_fce_bins)
    # Midpoints for plotting (optional)
    bin_centers = np.arange(len(hist_counts)) + 0.5  # centers at 0.5, 1.5, ..., N - 0.5
    # Plot with evenly spaced bars
    axsRight[1,1].bar(bin_centers, hist_counts, width=1.0, align='center')
    axsRight[1,1].set_xlabel(r'$f_{pe}/f_{ce}$')
    axsRight[1,1].set_xticks(np.arange(len(fpe_fce_bins)))
    axsRight[1,1].set_xticklabels([f"{val:.1f}" for val in fpe_fce_bins])

    hist_counts, _ = np.histogram(ae_vals[ae_vals<=1500.],bins=ae_bins)
    # Midpoints for plotting (optional)
    bin_centers = np.arange(len(hist_counts)) + 0.5  # centers at 0.5, 1.5, ..., N - 0.5
    axsRight[0,1].bar(bin_centers, hist_counts, width=1.0, align='center')
    axsRight[0,1].set_xlabel(r'$Ae\ index$')
    axsRight[0,1].set_xticks(np.arange(len(ae_bins)))
    axsRight[0,1].set_xticklabels([f"{val:.0f}" for val in ae_bins])


    axsLeft.set_title(r'$a)\ IQR_{norm}\ averages\ in\ Ae\ and\ f_{pe}/f_{ce}$', va='center', ha='center', fontsize=14)
    # Add side labels for each row
    axsRight[0,0].set_title(r'$b)\ Histogram\ of\ all\ Ae\ and\ f_{pe}/f_{ce}\ in\ linear\ bins$', va='center', ha='center', fontsize=14)
    
    axsRight[0,1].set_title(r'$c)\ Histogram\ of\ all\ Ae\ and\ f_{pe}/f_{ce}\ in\ quantile\ bins$', va='center', ha='center', fontsize=14)



    # now do the AE slices:
    # Get the tab20b colormap
    #cmap = plt.get_cmap('tab20b')
    # Get list of RGBA colors
    #colors = [cmap(i) for i in range(cmap.N)]

    #for i in range(10):

    #    axsSlice.plot(fpe_fce_centers,ae_fpe_array[i,:],label=f'Ae = {ae_bins[i]}-{ae_bins[i+1]}', color = colors[i])
    
    #axsSlice.set_xlabel(r'$f_{pe}/f_{ce}$')
    #axsSlice.set_ylabel(r'$IQR_{norm}$')
    #axsSlice.set_xticks(np.arange(len(fpe_fce_bins)))
    #axsSlice.set_xticklabels([f"{val:.2f}" for val in fpe_fce_bins])
    #axsSlice.set_xscale('log')
    #axsSlice.legend()

    plt.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig7/{name}_with_hists_low.png')
