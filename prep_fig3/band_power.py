# This code creates the data required to make the bands of power plots (low frequency, lowerband and upperband)

# CURRENTLY this is using slightly incomplete data as the input whilst the full job is running


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches

# other
import pandas as pd
import string

plt.style.use("ggplot")


def make_stat_histogram(dataframe,stat):
    
    # define spatial bins: 0.1 L and 1 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)

    power = stat#.tolist()
    dist_MLT_L = np.zeros(((len(abins)-1),(len(rbins)-1),len(power)))

    MLT = dataframe['MLT_rads'].values.tolist()
    L_star = dataframe['L_plot'].values.tolist()

    # Digitize the MLT and L values to get the bin indices
    MLT_bin_indices = np.digitize(MLT, abins) - 1  # Subtract 1 to make bins zero-indexed
    L_bin_indices = np.digitize(L_star, rbins) - 1

    # Create a 2D array to store the binned power data
    binned_power = np.zeros((len(abins)-1, len(rbins)-1))
    # create an array to count the np.nan's so that we can remove it from the toal count in a bin
    remove_bincounts = np.zeros((len(abins)-1, len(rbins)-1))
    # Aggregate power values by bin
    # note: want to track the distribution in each bin
    
    print(len(power))
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
    dist_MLT_L[dist_MLT_L==0.] = np.nan

    return binned_mean_power, dist_MLT_L


def make_dial_plot(ax,histogram,colorbar,title):  
  
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
chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/all_stats_mean_redidV12_fpefceADDED.nc')


# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"]

# Make a 4 x 3 plots (do: survey power, burst peak power, IQR, burst peak/survey)

# Create a figure

chorus_result = chorus_result.where((chorus_result["isChorus"]),drop=True)

# set conditions to just the equatorial - everything else!
fig, axes = plt.subplots(4, 3, subplot_kw={'projection': 'polar'}, figsize=(8, 10)) 
plot_names = list(string.ascii_lowercase)

for i,name in enumerate(names):

    print(f'{name} started...')
    conditions = (chorus_result[f'{name}_peak']<1e4)&(chorus_result["AE"]<9e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result["MLAT"]>=6.)

    survey_vals = np.where(conditions,chorus_result[f'{name}_survey_power'][:].values, np.nan)
    iqr_vals =  np.where(conditions,chorus_result[f'{name}_iqr'][:].values, np.nan)
    mad_vals = np.where(conditions,chorus_result[f'{name}_mad'][:].values, np.nan)
    peak_vals = np.where(conditions,chorus_result[f'{name}_peak'][:].values, np.nan)
    ratio = peak_vals/survey_vals
  
    
    # make histograms
    peak_hist,mlt_l_dists = make_stat_histogram(chorus_result,peak_vals)
    peakratio_hist, mlt_l_dists = make_stat_histogram(chorus_result,ratio)
    surv_hist, mlt_l_dists = make_stat_histogram(chorus_result,survey_vals)
    iqr_hist,mlt_l_dists = make_stat_histogram(chorus_result,iqr_vals)
    mad_hist,mlt_l_dists = make_stat_histogram(chorus_result,mad_vals)

    print(f'{name} histograms made...')
    # Kaine's smoothing technique before plotting

    survey_smooth = np.zeros((24*10,70))
    ratio_smooth = np.zeros((24*10,70))
    peak_smooth = np.zeros((24*10,70))
    iqr_smooth = np.zeros((24*10,70))
    mad_smooth= np.zeros((24*10,70))
    for k in range(24):
        survey_smooth[k*10:(k+1)*10,:] = surv_hist[k,:]

        ratio_smooth[k*10:(k+1)*10,:] = peakratio_hist[k,:]

        peak_smooth[k*10:(k+1)*10,:] = peak_hist[k,:]
        
        iqr_smooth[k*10:(k+1)*10,:] = iqr_hist[k,:]

        mad_smooth[k*10:(k+1)*10,:] = mad_hist[k,:]

    # save all data to list
    datas = [survey_smooth,peak_smooth,ratio_smooth,iqr_smooth]

    # define normalisations
    peak_norm = mcolors.LogNorm(vmin=10**(-4), vmax=10**(-2))
    survey_norm = mcolors.LogNorm(vmin=10**(-4), vmax=10**(-2))
    ratio_norm = mcolors.Normalize(vmin=1.,vmax=100.)
    iqr_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
    # save all to list
    norms = [survey_norm,peak_norm,ratio_norm,iqr_norm]

    cmaps = ['jet','jet','seismic','jet']

    # now do plot
    print(f'{name} plots begun...')
    
    for j,(data,norm) in enumerate(zip(datas,norms)):
        
        title = plot_names[(4*i)+j] +')'
        axis = axes[j,i]
        axis.set_title(title,loc='left',fontsize=14)

        if j == 2:
             dial_plot = make_dial_frac_plot(axis,data,norm,title)
        else:
            dial_plot = make_dial_plot(axis,data,norm,title)

        if (j == 0):
            if (i==0):
                axis.set_xticklabels([str(i) if (i==18 or i == 12) else '' for i in range(0,24)],fontsize=14)
            elif (i==1):
                axis.set_xticklabels([str(i) if (i == 12) else '' for i in range(0,24)],fontsize=14)
            else:
                axis.set_xticklabels([str(i) if (i == 12 or i==6) else '' for i in range(0,24)],fontsize=14)
        
        elif (j == 3):
            if (i==0):
                axis.set_xticklabels([str(i) if (i==18 or i == 0) else '' for i in range(0,24)],fontsize=14)
            elif (i==1):
                axis.set_xticklabels([str(i) if (i == 0) else '' for i in range(0,24)],fontsize=14)
            else:
                axis.set_xticklabels([str(i) if (i == 0 or i==6) else '' for i in range(0,24)],fontsize=14)
        
        else:

            if (i==0):
                axis.set_xticklabels([str(i) if (i==18) else '' for i in range(0,24)],fontsize=14)
            elif (i==2):
                axis.set_xticklabels([str(i) if (i==6) else '' for i in range(0,24)],fontsize=14)
            else:
                axis.set_xticklabels(['' for i in range(0,24)],fontsize=14)

        box = axis.get_position()
        axis.set_position([box.x0 - 0.1, box.y0, box.width, box.height])


    print(f'{name} done')


    # Add the color bar to the figure
    #cbar=fig.colorbar(sm, ax=axes[0,-1], orientation='vertical', fraction=0.1, pad=0.1)
    #cbar.ax.tick_params(labelsize=18)  # Set colorbar tick label font size
    #cbar.set_label(r'$nT$', fontsize=18) 
    #plt.subplots_adjust(hspace=0.1)
    #plt.tight_layout()
    #plt.gca().set_facecolor(color='#72A296')

    #plt.subplots_adjust(hspace=0.1,wspace = 0.01) 

# add the colorbar to the side

# first, for power:
sm = plt.cm.ScalarMappable(cmap='jet', norm=peak_norm)
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = fig.add_axes([0.85, 0.6, 0.05, 0.2])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar.set_label(r'$nT^2$', fontsize=14) 

# Add a separate axis for the color bar of the ratio
sm = plt.cm.ScalarMappable(cmap='Greens', norm=ratio_norm)
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax2 = fig.add_axes([0.85, 0.31, 0.05, 0.15])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax2.spines.values():
    spine.set_visible(False)
cbar_ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar2=fig.colorbar(sm, ax=cbar_ax2, orientation='vertical', fraction=0.5, pad=0.05,ticks=[10, 100])
cbar2.ax.tick_params(labelsize=14)  # Set colorbar tick label font size


# Add colorbar for the iqr
sm = plt.cm.ScalarMappable(cmap='jet', norm=iqr_norm)
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = fig.add_axes([0.85, 0.11, 0.05, 0.15])  # [left, bottom, width, height]
# Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
cbar=fig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar.set_label(r'$nT^2$', fontsize=14)


    
fig.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig3/HighLatV1.png')