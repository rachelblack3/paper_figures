
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
   

# Read the array from file
read_data = np.loadtxt('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/test.txt')
read_data = read_data.reshape((6,24*10,70))

survey_samp_smoothed = read_data[0,:,:]
burst_samp_smoothed = read_data[1,:,:]
burst_survey_sampling = read_data[2,:,:]
survey_chorus_occurrence = read_data[3,:,:]/read_data[0,:,:]
burst_chorus_occurrence = read_data[4,:,:]/read_data[1,:,:]
burst_survey_chorus = read_data[5,:,:]

# plot dial plot
print("starting plot...")
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 3, subplot_kw={'projection': 'polar'}, figsize=(10, 6)) 

title='a) Survey'
axis = axes[0,0]
axis.set_title(title,loc='center',fontsize=14)
sampling_norm = mcolors.LogNorm(vmin=10**(2), vmax=10**(5))
dial_plot = make_multiple_dial_plot(axis,survey_samp_smoothed,sampling_norm,title)
# first burst plot has 12 and 18 only

axis.set_xticklabels([str(i) if (i==18 or i == 12) else '' for i in range(0,24)],fontsize=14)
box1 = axes[0,0].get_position()
axes[0,0].set_position([box1.x0 - 0.1, box1.y0, box1.width, box1.height])

title = 'b) Burst'
axis = axes[0,1]
axis.set_title(title,loc='center',fontsize=14)
dial_plot = make_multiple_dial_plot(axis,burst_samp_smoothed,sampling_norm,title)
axis.set_xticklabels([str(i) if (i == 12) else '' for i in range(0,24)],fontsize=14)
# add in the mean power plot
# Manually move burst left, closer to survey
box2 = axes[0,1].get_position()
axes[0,1].set_position([box2.x0 - 0.14, box2.y0, box2.width, box2.height])

title = 'c) Burst/Survey'

ratio_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
axis = axes[0,2]
axis.set_title(title,loc='center',fontsize=14)
dial_plot = make_multiple_rel_dial_plot(axis,burst_survey_sampling,ratio_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==12) else '' for i in range(0,24)],fontsize=14)

axis = axes[1,0]
chorus_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
dial_plot = make_multiple_dial_plot(axis,survey_chorus_occurrence,chorus_norm,title)
# first burst plot has 12 and 18 only

axis.set_xticklabels([str(i) if (i==18 or i ==0) else '' for i in range(0,24)],fontsize=14)
box2 = axes[1,0].get_position()
axes[1,0].set_position([box2.x0 - 0.1, box2.y0, box2.width, box2.height])


axis = axes[1,1]
dial_plot = make_multiple_dial_plot(axis,burst_chorus_occurrence,chorus_norm,title)
axis.set_xticklabels([str(i) if (i ==0) else '' for i in range(0,24)],fontsize=14)
# add in the mean power plot
# Manually move burst left, closer to survey
box4 = axes[1,1].get_position()
axes[1,1].set_position([box4.x0 - 0.14, box4.y0, box4.width, box4.height])

ratio_norm = mcolors.LogNorm(vmin=10**(-2), vmax=10**(0))
axis = axes[1,2]
dial_plot = make_multiple_rel_dial_plot(axis,burst_survey_chorus,ratio_norm,title)
axis.set_xticklabels([str(i) if (i ==6 or i==0) else '' for i in range(0,24)],fontsize=14)
# Manually move burst left, closer to survey
box5 = axes[0,2].get_position()
axes[0,2].set_position([box5.x0 - 0.05, box5.y0, box5.width, box5.height])
box5 = axes[1,2].get_position()
axes[1,2].set_position([box5.x0 - 0.05, box5.y0, box5.width, box5.height])

# Add in the colorbars again

# Add a separate axis for the color bar of the SAMPLING
sm = plt.cm.ScalarMappable(cmap='jet', norm=mcolors.LogNorm(1e2,1e5))
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = fig.add_axes([0.5, 0.55, 0.05, 0.3])  # [left, bottom, width, height]
#Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
#
cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar.set_label(r'$N_{burst,survey}$', fontsize=14) 

# add seperate axis for colorbar of chorus occurrence
sm = plt.cm.ScalarMappable(cmap='jet', norm=mcolors.LogNorm(1e-2,1e0))
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = fig.add_axes([0.5, 0.15, 0.05, 0.3])  # [left, bottom, width, height]
#Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
#
cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar.set_label(r'$Occurrence_{burst,survey}$', fontsize=14)

# Add a separate axis for the color bar of the FRACTION

sm = plt.cm.ScalarMappable(cmap='Greens', norm=mcolors.LogNorm(1e-2,1e0))
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax2 = fig.add_axes([0.875, 0.37, 0.05, 0.3])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax2.spines.values():
    spine.set_visible(False)
cbar_ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar2=fig.colorbar(sm, ax=cbar_ax2, orientation='vertical', fraction=0.5, pad=0.05,ticks=[0.01,1,100])
cbar2.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar2.set_label(r'$N_{burst}/N_{survey}$', fontsize=14)



print("saving plot...")

fig.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/fig1V1.png')

