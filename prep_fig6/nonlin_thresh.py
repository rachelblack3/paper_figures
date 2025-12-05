# This code creates the data required to find the nonlinear threshold for all chorus events (low frequency, lowerband and upperband)

# CURRENTLY this is using slightly incomplete data as the input whilst the full job is running


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
# numerical essentials 
import numpy as np
import roman 

# for plotting
from matplotlib import pyplot as plt, animation
from metpy.plots import add_timestamp
import os
import pandas as pd

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches

# for cdf reading
from spacepy import toolbox
#spacepy.toolbox.update(leapsecs=True)
from spacepy import pycdf
import h5py

# other numerical tools
import os
import glob
from datetime import datetime,timedelta,date

import xarray as xr
import seaborn as sns
from scipy.stats import norm,iqr

plt.style.use("ggplot")



def make_stat_histogram(dataframe,stat):
    
    # define spatial bins: 0.1 L and 1 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)

    power = stat#.tolist()
    dist_MLT_L = np.zeros(((len(abins)-1),(len(rbins)-1),len(power)))

    MLT = dataframe['MLT_rads'].values.tolist()
    L_star = dataframe['Lstar'].values.tolist()

    # Digitize the MLT and L values to get the bin indices
    MLT_bin_indices = np.digitize(MLT, abins) - 1  # Subtract 1 to make bins zero-indexed
    L_bin_indices = np.digitize(L_star, rbins) - 1

    # Create a 2D array to store the binned power data
    binned_power = np.zeros((len(abins)-1, len(rbins)-1))
    # create an array to count the np.nan's so that we can remove it from the toal count in a bin
    remove_bincounts = np.zeros((len(abins)-1, len(rbins)-1))
    # Aggregate power values by bin
    # note: want to track the distribution in each bin
    

    for i in range(len(power)):
        mlt_idx = MLT_bin_indices[i]
        l_idx = L_bin_indices[i]
        if 0 <= mlt_idx < len(abins)-1 and 0 <= l_idx < len(rbins)-1:  # Ensure valid indices


            if np.isnan(power[i]):
                remove_bincounts[mlt_idx, l_idx] +=1
                dist_MLT_L[mlt_idx,l_idx,i] = np.nan
            elif power[i]>1*10**3:
                remove_bincounts[mlt_idx, l_idx] +=1
                dist_MLT_L[mlt_idx,l_idx,i] = np.nan
            else:
                binned_power[mlt_idx, l_idx] += power[i]
                dist_MLT_L[mlt_idx,l_idx,i]=power[i]
            
   
    # Optionally, calculate the mean by counting occurrences in each bin
    bin_counts = np.histogram2d(MLT, L_star, bins=[abins, rbins])[0]
    bin_counts = bin_counts-remove_bincounts
    
    
    #for count in bin_counts:
        #print(count)
    binned_mean_power = np.divide(binned_power, bin_counts, where=bin_counts > 10)  # Avoid division by zero


    # To replace the result with NaN where the condition isn't met:
    binned_mean_power = np.where(bin_counts > 10, binned_mean_power, np.nan)
   
    # Print the inf mask
    # Replace infinities with a specific value (e.g., 0)
    
    # 'binned_mean_power' now holds the average power in each MLT-L bin
    #dist_MLT_L[dist_MLT_L==0.] = np.nan

    all_bin_counts = np.sum(bin_counts)
   
    return binned_mean_power, dist_MLT_L,all_bin_counts


def make_dial_plot(ax,histogram,colorbar,title,cmap):  
  
# Plot dial plot

    L_list = np.linspace(0, 7, 1)
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 70+1,endpoint=True)
    abins = np.linspace(0,2*np.pi, (24*10)+1,endpoint=True)
    A, R = np.meshgrid(abins, rbins)

    #ax.set_rticks(L_list)  # Less radial ticks
    #ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    #Set the color of radial labels
    # Set custom y-tick positions and labels
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
    ax.set_yticklabels(['', '', '', '2', '3', '', '','7'])

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2 * np.pi, 24, endpoint=False))
    #ax.set_xticklabels([])
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)],fontsize=8)
    funct= np.transpose(histogram)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#ffffff","#00dca5","#00AB7E","#022020"])
    colorbar_norm = colorbar
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap=cmap,edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    
    
    # Add radial dotted lines
    ax.grid(which='both', color='darkgrey', linestyle='-', linewidth=1, alpha=1, zorder=1) 

    # Add half black half white circle to the center
    circle_radius = 1 # Adjust radius as needed

    #ax.set_title(f"{title}",fontsize=6)
    #plt.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/plots/2013_2019/AE/out_high.png')
    # Adjust layout to prevent overlap
    # Create semicircles
    semicircle_white = patches.Wedge((0, 0), circle_radius, 0, 180, facecolor='white', edgecolor='black', transform=ax.transData._b, zorder=10)
    semicircle_black = patches.Wedge((0, 0), circle_radius, 180, 360, facecolor='black', edgecolor='black', transform=ax.transData._b,zorder=10)

    # Add semicircles to the plot
    ax.add_patch(semicircle_white)
    ax.add_patch(semicircle_black)

    # Add semicircles to the plot
    ax.add_patch(semicircle_white)
    ax.add_patch(semicircle_black)
    #plt.figure(figsize=(10,10))



def make_dial_frac_plot(ax,histogram,colorbar,title):  
  
# Plot dial plot

    L_list = np.linspace(0, 7, 1)
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 70+1,endpoint=True)
    abins = np.linspace(0,2*np.pi, (24*10)+1,endpoint=True)
    A, R = np.meshgrid(abins, rbins)

    #ax.set_rticks(L_list)  # Less radial ticks
    #ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    #Set the color of radial labels
    # Set custom y-tick positions and labels
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
    ax.set_yticklabels(['', '', '', '2', '3', '', '','7'])

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2 * np.pi, 24, endpoint=False))
    #ax.set_xticklabels([])
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)],fontsize=8)
    funct= np.transpose(histogram)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#ffffff","#00dca5","#00AB7E","#022020"])
    colorbar_norm = colorbar
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap='Greens',edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    
    
    # Add radial dotted lines
    ax.grid(which='both', color='darkgrey', linestyle='-', linewidth=1, alpha=1, zorder=1) 

    # Add half black half white circle to the center
    circle_radius = 1 # Adjust radius as needed

    #ax.set_title(f"{title}",fontsize=6)
    #plt.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/plots/2013_2019/AE/out_high.png')
    # Adjust layout to prevent overlap
    # Create semicircles
    semicircle_white = patches.Wedge((0, 0), circle_radius, 0, 180, facecolor='white', edgecolor='black', transform=ax.transData._b, zorder=10)
    semicircle_black = patches.Wedge((0, 0), circle_radius, 180, 360, facecolor='black', edgecolor='black', transform=ax.transData._b,zorder=10)

    # Add semicircles to the plot
    ax.add_patch(semicircle_white)
    ax.add_patch(semicircle_black)

    # Add semicircles to the plot
    ax.add_patch(semicircle_white)
    ax.add_patch(semicircle_black)
    #plt.figure(figsize=(10,10))


# import the chorus dataset
chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')

# name the different bands for loop
names = ["lf", "lower", "upper"]
central_freq = [0.05, 0.3, 0.7]
labels = [roman.toRoman(i).lower() for i in range(1, 19)]

chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

# Make a 3 x 1 plots (do: non/quasi linear for each band)
# Create a figure
fig,axes = plt.subplots(1, 3, subplot_kw={'projection': 'polar'}, figsize=(12, 6))
fig = plt.figure(layout='constrained', figsize=(20, 10))
subfigs = fig.subfigures(2, 1)
axsBurstThresh = subfigs[0].subplots(1, 3, subplot_kw={'projection': 'polar'}) 

for i,name in enumerate(names):

    conditions = (chorus_result["Lstar"]>2.)
    non_quas_vals = np.where(conditions,chorus_result[f'thresh_{name}20'][:].values,np.nan)
    # set conditions
    non_quas,mlt_l_dists,all_bin_counts = make_stat_histogram(chorus_result,non_quas_vals)
    nonlinear_ratio = np.zeros((24*10,70))
   
    for k in range(24):
        nonlinear_ratio[k*10:(k+1)*10,:] = non_quas[k,:]

    ratio_norm = mcolors.Normalize(vmin=0., vmax=0.1)

    print(np.nanmax(nonlinear_ratio),np.nanmin(nonlinear_ratio))
    dial_plot = make_dial_frac_plot(axsBurstThresh[i],nonlinear_ratio,ratio_norm,f'{name}-band nonlinear ratio')

axsBurstThresh[0].set_xticklabels([str(i) if (i==12 or i==18) else '' for i in range(0,24)],fontsize=14)
axsBurstThresh[1].set_xticklabels([str(i) if (i==12) else '' for i in range(0,24)],fontsize=14)
axsBurstThresh[2].set_xticklabels([str(i) if (i==12 or i==6) else '' for i in range(0,24)],fontsize=14)

for ax in axsBurstThresh:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])


# Add colorbar
sm = plt.cm.ScalarMappable(cmap='Greens', norm=ratio_norm)
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = subfigs[0].add_axes([0.9, 0.3, 0.01, 0.4])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar=subfigs[0].colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.8, pad=0.05,ticks=[0.,0.1])
cbar.ax.tick_params(labelsize=14)


axsSurvThresh = subfigs[1].subplots(1, 3, subplot_kw={'projection': 'polar'}) 

for i,name in enumerate(names):

    conditions = (chorus_result["Lstar"]>2.)
    non_quas_vals = np.where(conditions,chorus_result[f'thresh_{name}_survey_20'][:].values,np.nan)
    # set conditions
    non_quas,mlt_l_dists,all_bin_counts = make_stat_histogram(chorus_result,non_quas_vals)

    nonlinear_ratio = np.zeros((24*10,70))
   
    for k in range(24):
        nonlinear_ratio[k*10:(k+1)*10,:] = non_quas[k,:]

    ratio_norm = mcolors.Normalize(vmin=0., vmax=0.1)

    print(np.nanmax(nonlinear_ratio),np.nanmin(nonlinear_ratio))
    dial_plot = make_dial_frac_plot(axsSurvThresh[i],nonlinear_ratio,ratio_norm,f'{name}-band nonlinear ratio')

axsSurvThresh[0].set_xticklabels([str(i) if (i==0 or i==18) else '' for i in range(0,24)],fontsize=14)
axsSurvThresh[1].set_xticklabels([str(i) if (i==0) else '' for i in range(0,24)],fontsize=14)
axsSurvThresh[2].set_xticklabels([str(i) if (i==0 or i==6) else '' for i in range(0,24)],fontsize=14)

for ax in axsSurvThresh:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Add colorbar
sm = plt.cm.ScalarMappable(cmap='Greens', norm=ratio_norm)
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = subfigs[1].add_axes([0.9, 0.3, 0.01, 0.4])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar=subfigs[1].colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.8, pad=0.05,ticks=[0.,0.1])
cbar.ax.tick_params(labelsize=14)

fig.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig6/20deg_surveyComp010725.png')