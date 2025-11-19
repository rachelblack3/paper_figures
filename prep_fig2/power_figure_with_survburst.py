import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import xarray as xr
import pandas as pd
# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches

# functions for making the power dial plots
plt.style.use("ggplot")

def make_histogram(dataframe):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)
    A, R = np.meshgrid(abins, rbins)


    # find BURST histogram first
    mlt_l_array = np.array(list(zip(dataframe['MLT_rads'].values.tolist(), dataframe['Lstar'].values.tolist())))
    #calculate histogram
    hist,_,_ = np.histogram2d(mlt_l_array[:,0],mlt_l_array[:,1],bins=(abins, rbins))

    return hist

def make_stat_histogram(dataframe,stat):
    
    # define spatial bins: 0.1 L and 1 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)

    power = stat.tolist()
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
    binned_mean_power = np.divide(binned_power, bin_counts, where=bin_counts > 50)  # Avoid division by zero
    # To replace the result with NaN where the condition isn't met:
    binned_mean_power = np.where(bin_counts > 50, binned_mean_power, np.nan)
   
    # Print the inf mask
    # Replace infinities with a specific value (e.g., 0)
    
    # 'binned_mean_power' now holds the average power in each MLT-L bin
    dist_MLT_L[dist_MLT_L==0.] = np.nan

    return binned_mean_power, dist_MLT_L


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
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap='plasma',edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    
    
    # Add radial dotted lines
    ax.grid(which='both', color='darkgrey', linestyle='-', linewidth=1, alpha=1, zorder=1) 

    # Add half black half white circle to the center
    circle_radius = 1 # Adjust radius as needed


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



# also read in the data from the last figure on sampling: low MLAT and high MLAT
#
names = ["Survey_sampling","Burst_sampling","Burst_survey_sampling","Survey_chorus","Burst_chorus","Burst_survey_chorus"]
arrays = [np.loadtxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/VanAllenB/LowMLATdataV1_{name}.txt') for name in names]
read_data_occurrence = np.stack(arrays)  # shape (6,240,70)

# read survey and burst chorus occurrence
survey_chorus_occurrence = read_data_occurrence[3,:,:]/read_data_occurrence[0,:,:]
burst_chorus_occurrence = read_data_occurrence[4,:,:]/read_data_occurrence[1,:,:]

names = ["Survey_power_just_chorus_events","Survey_power_all_events","burst_power_just_chorus_events","burst_power_all_events","survey_for_burst"]
arrays = [np.loadtxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig2/VanAllenB/LowMLATdataV1_{name}.txt') for name in names]
read_data = np.stack(arrays)  # shape (6,240,70)
#names = ["Survey_sampling","Burst_sampling","Burst_survey_sampling","Survey_chorus","Burst_chorus","Burst_survey_chorus"]
#arrays = [np.loadtxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/VanAllenB/LowMLATdataV1_{name}.txt') for name in names]
#read_data = np.stack(arrays)  # shape (6,240,70)

# Multiply the power by the chorus occurrence, so that we can get the effective power (?maybe different name)
survey_power_smooth = read_data[0,:,:]
relative_survey_power = read_data[1,:,:]
burst_power_smooth = read_data[2,:,:]
relative_burst_power = read_data[3,:,:]
survey_for_burst = read_data[4,:,:]
relative_survey_for_burst = survey_for_burst*burst_chorus_occurrence

chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/burstB270725_justChorus.nc')#xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus_burst260525.nc')
# save total power to chorus_result
burst_pow = chorus_result["mean_burst_int"].where((np.abs(chorus_result["MLAT"])<6.))
print("burst read...")

survey = survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerB_Version20250827.csv')#pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerA_Version20240515_withchorus_pos.csv')
# find the chorus in survey 
survey_chorus = survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>2.) & (np.abs(survey["MLAT"])<6.)]
survey_pow = survey_chorus["Survey_integral"]

print("plotting LHS...")

fig = plt.figure(layout='constrained', figsize=(18, 5.33))
subfigs = fig.subfigures(1, 2, wspace=0.07)

# do plot

axsLeft = subfigs[0].subplots(2, 3, subplot_kw={'projection': 'polar'}) 
title=r'$a)$'
axis = axsLeft[0,0]
axis.set_title(title,loc='Left',fontsize=14)
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,survey_power_smooth,power_norm,title)
axis.set_xticklabels([str(i) if (i==18 or i == 12) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0 - 0.01, box.y0, box.width, box.height])

# burst power (just chor events)
axis = axsLeft[0,1]
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,burst_power_smooth,power_norm,title)
axis.set_xticklabels([str(i) if (i == 12) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0, box.y0, box.width, box.height])

axis = axsLeft[0,2]
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,survey_for_burst,power_norm,title)
axis.set_xticklabels([str(i) if (i==6 or i == 12) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0 - 0.01, box.y0, box.width, box.height])

title=r'$b)$'
axis = axsLeft[1,0]
axis.set_title(title,loc='Left',fontsize=14)
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,relative_survey_power,power_norm,title)
axis.set_xticklabels([str(i) if (i==18 or i==0) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0 - 0.01, box.y0, box.width, box.height])


# burst power (all events - so x occurrence)
axis = axsLeft[1,1]
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,relative_burst_power,power_norm,title)
axis.set_xticklabels([str(i) if (i==0) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0, box.y0, box.width, box.height])


# burst power (just chor events)
axis = axsLeft[1,2]
power_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
dial_plot = make_multiple_dial_plot(axis,relative_survey_for_burst,power_norm,title)
axis.set_xticklabels([str(i) if (i==6 or i==0) else '' for i in range(0,24)],fontsize=14)
box = axis.get_position()
axis.set_position([box.x0, box.y0, box.width, box.height])

sfig = subfigs[0]
#Add a separate axis for the color bar of the SAMPLING
sm = plt.cm.ScalarMappable(cmap='plasma', norm=mcolors.LogNorm(1e-6,1e-2))
sm.set_array([])  # Required for the colorbar to work properly
cbar_ax = sfig.add_axes([0.92, 0.3, 0.05, 0.4])  # [left, bottom, width, height]
cbar_ax.set_facecolor('none')
#Make the axes invisible
for spine in cbar_ax.spines.values():
    spine.set_visible(False)
cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# Add the color bar to the figure
#
cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
cbar.set_label(r'$Power\ (nT^2)$', fontsize=14) 

axsLeft[0,0].text(0.5, 1.25, "i. survey", transform=axsLeft[0,0].transAxes,
        ha='center', fontsize=12)

axsLeft[0,1].text(0.5, 1.25, "ii. burst", transform=axsLeft[0,1].transAxes,
        ha='center', fontsize=12)

axsLeft[0,2].text(1.5, 1.25, "iii. survey for burst subset", transform=axsLeft[0,1].transAxes,
        ha='center', fontsize=12)

# Now do the right hand side 
print("plotting RHS...")


axsRight = subfigs[1].subplots(1, 1,) 
title=r'$c)$'
# survey rel
log_bins = np.linspace(-8, 0, 80)
#sns.histplot(data=relative_survey_power.flatten(),bins=log_bins,color='darkorange',label="Survey", kde=False,stat='probability')

# burst rel
#log_bins = np.logspace(np.log10(np.nanmin(relative_burst_power)), np.log10(np.nanmax(relative_burst_power)), 50)

#sns.histplot(np.log10(survey_power_smooth.flatten()),stat='probability',binwidth=0.1,color='darkgreen', label="Survey")
#sns.histplot(np.log10(burst_power_smooth.flatten()),stat='probability',binwidth=0.1,color='lightgreen', label="Burst")

sns.histplot(np.log10(survey_pow),stat='probability',bins=log_bins,color='darkgreen', label="Survey")
sns.histplot(np.log10(burst_pow),stat='probability',bins=log_bins,color='lightgreen', label="Burst")

#plt.hist(all_dists,bins=log_bins)
# Add KDE for log-transformed data
axsRight.set_title(title,loc='left',fontsize=14)
#axsRight.set_xscale('log')
sns.set(font_scale=0.8)
axsRight.set_ylabel(r'$Probabillity$',fontsize=14)
axsRight.set_xlabel(r'$Chorus\ average\ power\ (nT^2)$',fontsize=14)
axsRight.set_xlim(-8,0)
axsRight.legend()
#plt.style.use('ggplot')
#plt.figure(figsize=(6,5))
#plt.title(r'$Integrated\ chorus\ power\ (f_{lhr} - f_{ce})$',fontsize=14)
#burst_pows = np.array(burst_chor["total_power"]).flatten()
#sns.histplot(np.log10(survey_chorus["Survey_integral"]),stat='probability',binwidth=0.1,color='darkgreen',label='Survey chorus events')
#burst_pows = burst_pows[burst_pows<10**5]
#sns.histplot(np.log10(burst_pows),stat='probability',binwidth=0.1,color='lightgreen', label="Burst chorus events")
#plt.xlabel(r'$Power\ (nT^{2})$',fontsize=14)
#axsRight.set_xlim(-6,0)
##plt.ylabel(r'$Probability$',fontsize=14)
#plt.legend()

# Shift row 2 (bottom row) to the right
for j in range(2):
    box = axsLeft[j, 0].get_position()
    axsLeft[j, 0].set_position([box.x0 + 0.02, box.y0, box.width, box.height])

# Add side labels for each row
sfig.text(0.04, 0.71, 'Average over\nchorus events', va='center', ha='center', fontsize=14, rotation='vertical')
sfig.text(0.04, 0.27, 'Average over\nall events', va='center', ha='center', fontsize=14, rotation='vertical')


plt.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig2/vanallenB/LowattemptB.png')
print("Saved...")
