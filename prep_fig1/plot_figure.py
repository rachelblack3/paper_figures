
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.ticker import ScalarFormatter
from datetime import datetime

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches
import xarray as xr
from datetime import timedelta


def make_multiple_dial_plot(ax,histogram,colorbar,title):  
  
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
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap='jet',edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    
    
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

def make_multiple_rel_dial_plot(ax,histogram,colorbar,title):  
  
# Plot dial plot

    L_list = np.linspace(0, 7, 1)
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 70+1,endpoint=True)
    abins = np.linspace(0,2*np.pi, (24*10)+1,endpoint=True)
    A, R = np.meshgrid(abins, rbins)


    # plot it!
    #ax.set_rticks(L_list)  # Less radial ticks
    #ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    #Set the color of radial labels
    for label in ax.get_yticklabels():
        label.set_color('black')  # Change 'red' to any color you prefer
        label.set_fontweight('bold')

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2 * np.pi, 24, endpoint=False))
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)])
    funct= np.transpose(histogram)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#ffffff","#00dca5","#00AB7E","#022020"])
    colorbar_norm = colorbar
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap='Greens',edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    #cbar=fig.colorbar(pc,label=r'$N_{Burst}/N_{Survey}$',pad=0.1)

    # Add radial dotted lines
    ax.grid(which='both', color='darkgrey', linestyle='-', linewidth=1, alpha=1, zorder=1) 

    # Add half black half white circle to the center
    circle_radius = 1 # Adjust radius as needed

   
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
   

# Read the array from disk
read_data = np.loadtxt('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/test.txt')
read_data = read_data.reshape((6,24*10,70))

survey_samp_smoothed = read_data[0,:,:]
burst_samp_smoothed = read_data[1,:,:]
burst_survey_sampling = read_data[2,:,:]
survey_chorus_occurrence = read_data[3,:,:]
burst_chorus_occurrence = read_data[4,:,:]
burst_survey_chorus = read_data[5,:,:]

# plot dial plot
print("starting plot...")
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 3, subplot_kw={'projection': 'polar'}, figsize=(12, 4)) 
title = 'Fraction of burst-mode triggers to total survey time'

title='a) survey sampling'
axis = axes[0,0]
axis.set_title(title,loc='left',fontsize=14)
sampling_norm = mcolors.LogNorm(vmin=10**(2), vmax=10**(5))
dial_plot = make_multiple_dial_plot(axis,survey_samp_smoothed,sampling_norm,title)
# first burst plot has 12 and 18 only

axis.set_xticklabels([str(i) if (i==18 or i == 12 or i ==6) else '' for i in range(0,24)],fontsize=14)

title = 'b) burst sampling'

axis = axes[0,1]
axis.set_title(title,loc='left',fontsize=14)
dial_plot = make_multiple_dial_plot(axis,burst_samp_smoothed,sampling_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==18) else '' for i in range(0,24)],fontsize=14)
# add in the mean power plot


title = 'c) burst/survey sampling'

ratio_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
axis = axes[0,2]
axis.set_title(title,loc='left',fontsize=14)
dial_plot = make_multiple_rel_dial_plot(axis,burst_survey_sampling,ratio_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==18 or i==0) else '' for i in range(0,24)],fontsize=14)


title='d) chorus occurrence in survey'
axis = axes[1,0]
axis.set_title(title,loc='left',fontsize=14)
chorus_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
dial_plot = make_multiple_dial_plot(axis,survey_chorus_occurrence,chorus_norm,title)
# first burst plot has 12 and 18 only

axis.set_xticklabels([str(i) if (i==18 or i == 12 or i ==6) else '' for i in range(0,24)],fontsize=14)

title = 'e) chorus occurrence in burst'

axis = axes[1,1]
axis.set_title(title,loc='left',fontsize=14)
dial_plot = make_multiple_dial_plot(axis,burst_chorus_occurrence,chorus_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==18) else '' for i in range(0,24)],fontsize=14)
# add in the mean power plot


title = 'f) burst/survey chorus occurrence'

ratio_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
axis = axes[1,2]
axis.set_title(title,loc='left',fontsize=14)
dial_plot = make_multiple_rel_dial_plot(axis,burst_survey_chorus,ratio_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==18 or i==0) else '' for i in range(0,24)],fontsize=14)



fig.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/fig1V1.png')

