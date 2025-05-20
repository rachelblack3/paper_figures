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


# functions

def to_dt(nanoseconds_value):

    # Convert nanoseconds to seconds
    seconds_value = nanoseconds_value / 1e9

    # Unix epoch
    epoch = datetime(1970, 1, 1)

    # Convert to datetime
    human_readable_datetime = epoch + timedelta(seconds=seconds_value)

    return human_readable_datetime

# Flatten lists function

def flattenColumn(input, column):
    '''
    column is a string of the column's name.
    for each value of the column's element (which might be a list),
    duplicate the rest of columns at the corresponding row with the (each) value.
    '''
    column_flat = pd.DataFrame(
        [
            [i, c_flattened]
            for i, y in input[column].apply(list).iteritems()
            for c_flattened in y
        ],
        columns=['I', column]
    )
    column_flat = column_flat.set_index('I')
    return (
        input.drop(column, 1)
             .merge(column_flat, left_index=True, right_index=True)
    )


# Function to convert string representation of list to actual list of floats
def string_to_float_list(s):
    # Remove the brackets and split the string by spaces
    s = s.strip('[]')
    # Split by spaces, filter out empty strings, and convert to floats
    return [float(x) for x in s.split() if x]


def make_histogram_da(dataframe):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)
    A, R = np.meshgrid(abins, rbins)


    # find BURST histogram first
    mlt_l_array = np.array(list(zip(dataframe['MLT_rads'], dataframe['Lstar'])))
    #calculate histogram
    hist,_,_ = np.histogram2d(mlt_l_array[:,0],mlt_l_array[:,1],bins=(abins, rbins))

    return hist

def make_histogram(dataframe):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)
    A, R = np.meshgrid(abins, rbins)


    # find BURST histogram first
    mlt_l_array = np.array(list(zip(dataframe['MLT_rads'].tolist(), dataframe['Lstar'].tolist())))
    #calculate histogram
    hist,_,_ = np.histogram2d(mlt_l_array[:,0],mlt_l_array[:,1],bins=(abins, rbins))

    return hist


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
    

# starting with survey, access file that has all survey sampling + chorus occurrence

survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerA_Version20240515_withchorus_pos.csv')
print("survey loaded...")

# change all of the timestamps to datetime objects
survey["Timestamp"] = survey["Timestamp"].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
# change MLT to radians for binning in polar coords
survey["MLT_rads"] = 2*np.pi*survey["MLT"]/24
print("survey time and MLT loaded...")

# create AE and plasmapause filters
survey_low = survey[(survey["AE"]<100)]
survey_medium = survey[(survey["AE"]>=100)&(survey["AE"]<300)]
survey_high =survey[(survey["AE"]>=300)]

# keep only between 2013-2019 for now
time_condition = (survey["Timestamp"] > datetime(year=2012,month=12,day=31)) & (survey["Timestamp"] < datetime(year=2019,month=1,day=1))
survey = survey[time_condition]

# SAMPLING
# make histogram
survey_samp_hist = make_histogram(survey)
# remove zeroes
zero_mask = survey_samp_hist == 0
survey_samp_hist[zero_mask]=np.nan
# Kaine's smoothing technique
survey_samp_smoothed = np.zeros((24*10,70))
for i in range(24):
    survey_samp_smoothed[i*10:(i+1)*10,:] = survey_samp_hist[i,:]

# also remove bins for which there was low sampling for survey, which would skew the statistics. We have chosen 1000 samples as lowest
samp_mask =  (survey_samp_smoothed < 0.1*np.nanmedian(survey_samp_smoothed))
# Print the inf mask
print("\nLow sampling bin mask:")
print(samp_mask)
survey_samp_smoothed[samp_mask] = np.nan


# CHORUS EVENTS
# make histogram of chorus
survey_chorus_hist = make_histogram(survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>2.)])
# remove zeroes
zero_mask = survey_chorus_hist == 0
survey_chorus_hist[zero_mask]=np.nan
# Kaine's smoothing technique
survey_chorus_smoothed = np.zeros((24*10,70))
for i in range(24):
    survey_chorus_smoothed[i*10:(i+1)*10,:] = survey_chorus_hist[i,:]

# also remove bins for which there was low sampling for survey, which would skew the statistics. We have chosen 1000 samples as lowest
samp_mask_survey =  (survey_chorus_smoothed < 0.1*np.nanmedian(survey_chorus_smoothed))
# Print the inf mask
print("\nLow sampling bin mask:")
print(samp_mask_survey)
survey_chorus_smoothed[samp_mask_survey] = np.nan

print("survey dists made...")


# CHORUS OCCURRENCE (chorus events/all events)
survey_chorus_occurrence = survey_chorus_smoothed/survey_samp_smoothed



# Now doing burst
burst = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/curr_combined/full_13_18.nc')
# Reduce 'chorus_flag' over the 'y' dimension, grouped by 'x'
is_chorus = burst['chorus_flag'].any(dim='y')  # Check for any True in 'y'
burst['isChorus'] = is_chorus
print("burst loaded...")

# create AE and plasmapause filters
#burst_low = burst.where(burst["AE"]<100)
#burst_medium= burst.where((burst["AE"]>=100)&(burst["AE"]<300))
#data_high =burst.where(burst["AE"]>=300)

# change MLT to radians for binning in polar coords
burst["MLT_rads"] = 2*np.pi*burst["MLT"]/24

# SAMPLING
# make histogram
burst_samp_hist = make_histogram_da(burst)
print("burst sampling histogram made...")
# remove zeroes
zero_mask = burst_samp_hist == 0
burst_samp_hist[zero_mask]=np.nan

# Kaine's smoothing technique
burst_samp_smoothed = np.zeros((24*10,70))
for i in range(24):
    burst_samp_smoothed[i*10:(i+1)*10,:] = burst_samp_hist[i,:]


# CHORUS EVENTS
# make histogram
burst_chorus_hist = make_histogram_da(burst.where((burst["isChorus"]==True) & (burst["Lstar"]>2.)))
print("burst chorus histogram made...")

# remove zeroes
zero_mask = burst_chorus_hist == 0
burst_chorus_hist[zero_mask]=np.nan

# Kaine's smoothing technique
burst_chorus_smoothed = np.zeros((24*10,70))
for i in range(24):
    burst_chorus_smoothed[i*10:(i+1)*10,:] = burst_chorus_hist[i,:]

burst_chorus_smoothed[samp_mask_survey] = np.nan

# CHORUS OCCURRENCE (chorus events/all events)
burst_chorus_occurrence = burst_chorus_smoothed/burst_samp_smoothed
print("burst dists made...")

# BURST/SURVEY COMPARISONS
#1. SAMPLING
burst_survey_sampling = burst_samp_smoothed/survey_samp_smoothed

#2. CHORUS OCCURRENCE
burst_survey_chorus = burst_chorus_smoothed/survey_chorus_smoothed

# plot all 6

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

