# Code that sums up power in the 3 'population' frequency bands, saves to:
# exohiss_power, lowerband_power, upper_bandpower variables

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
# numerical essentials 
import numpy as np

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
from scipy.stats import norm,iqr,median_abs_deviation

def to_dt(nanoseconds_value):

    # Convert nanoseconds to seconds
    seconds_value = nanoseconds_value / 1e9

    # Unix epoch
    epoch = datetime(1970, 1, 1)

    # Convert to datetime
    human_readable_datetime = epoch + timedelta(seconds=seconds_value)

    return human_readable_datetime


# functions
def skewness(data):
    N = len(data)
    std = np.std(data)
    mean = np.mean(data)

    skew_coeff = np.sqrt(N*(N-1))/(N-2)

    skew_coeff= skew_coeff/std**3

    skew=0
    for i in range(N):

        skew = skew + (data[i]-mean)**3/N
    
    return skew_coeff*skew

def second_coeff(data):

    N = len(data)
    std = np.nanstd(data)
    mean = np.nanmean(data)
    median = np.nanmedian(data)

    skew_coeff = 3

    skew = (mean-median)/std
    
    return skew_coeff*skew



def make_histogram(dataframe):
    
    # define spatial bins: 0.1 L and 0.25 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)
    A, R = np.meshgrid(abins, rbins)


    # find BURST histogram first
    mlt_l_array = np.array(list(zip(dataframe['MLT_rads'].values.tolist(), dataframe['L_plot'].values.tolist())))
    #calculate histogram
    hist,_,_ = np.histogram2d(mlt_l_array[:,0],mlt_l_array[:,1],bins=(abins, rbins))

    return hist

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
            elif power[i]==0:
                remove_bincounts[mlt_idx, l_idx] +=1
                dist_MLT_L[mlt_idx,l_idx,i] = np.nan
            else:
                binned_power[mlt_idx, l_idx] += power[i]
                dist_MLT_L[mlt_idx,l_idx,i]=power[i]

    # Optionally, calculate the mean by counting occurrences in each bin
    bin_counts = np.histogram2d(MLT, L_star, bins=[abins, rbins])[0]
    #bin_counts = bin_counts-remove_bincounts
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
    ax.set_yticklabels(['', '', '1', '', '', '', '','7'])

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2 * np.pi, 24, endpoint=False))
    #ax.set_xticklabels([])
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)],fontsize=18)
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



burst = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/all_burst230525.nc')
# Reduce 'chorus_flag' over the 'y' dimension, grouped by 'x'
is_chorus = burst['chorus_flag'].any(dim='y')  # Check for any True in 'y'
burst['isChorus'] = is_chorus

chorus_result = burst.where((burst["isChorus"]),drop=True)

# name the variables for holding total_power in different bands - survey
exohiss_survey_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="exohiss_survey_power"
)
lowerband_survey_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="lowerband_survey_power"
)
upperband_survey_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="upperband_survey_power"
)

# name the variables for holding total_power in different bands - burst
exohiss_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="exohiss_power"
)
lowerband_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="lowerband_power"
)
upperband_power = xr.DataArray(
    data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="upperband_power"
)

# name the variables for holding power peaks in different bands
exohiss_peak = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
        dims="x",
        name="exohiss_peak"
    )

lowerband_peak = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
        dims="x",
        name="lowerband_pek"
    )
upperband_peak = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 2D: (x, z)
        dims="x",
        name="upperband_peak"
    )

# name the variables for holding power means in different bands
exohiss_mean = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="exohiss_mean"
    )

lowerband_mean = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="lowerband_mean"
    )

upperband_mean = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="upperband_mean"
    )

# name the variables for holding power medians in different bands
exohiss_median= xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="exohiss_median"
    )

lowerband_median = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="lowerband_median"
    )

upperband_median = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="upperband_median"
    )

# name the variables for holding power means in different bands
exohiss_range = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="exohiss_range"
    )

lowerband_range = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="lowerband_range"
    )

upperband_range = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name="upperband_range"
    )

exohiss_iqr = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'exohiss_IQR'
    )

lowerband_iqr = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'lowerband_IQR'
    )

upperband_iqr = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'upperband_IQR'
    )

exohiss_mad = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'exohiss_mad'
    )

lowerband_mad = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'lowerband_mad'
    )

upperband_mad = xr.DataArray(
        data=np.full(burst.sizes["x"], np.nan),  # 1D: (x)
        dims="x",
        name=f'upperband_mad'
    )

band_ranges = [[0,1],[1,5],[5,-1]]
bands_survey = [exohiss_survey_power,lowerband_survey_power,upperband_survey_power]
bands_burst = [exohiss_power,lowerband_power,upperband_power]

peak_bands = [exohiss_peak,lowerband_peak,upperband_peak]
mean_bands = [exohiss_mean,lowerband_mean,upperband_mean]
median_bands = [exohiss_median,lowerband_median,upperband_median]
range_bands = [exohiss_range,lowerband_range,upperband_range]
iqr_bands = [exohiss_iqr,lowerband_iqr,upperband_iqr]
mad_bands = [exohiss_mad,lowerband_mad,upperband_mad]
names = ["exohiss_band", "lowerband", "upperband"]

for variable_survey, variable_burst, variable_peak, variable_IQR, variable_mad, variable_mean, variable_median, band_range, name in zip(bands_survey,bands_burst,peak_bands,iqr_bands, mad_bands, mean_bands, median_bands, band_ranges, names):
    print(f'starting {name}...')
    print(f'Number of samples: {burst.sizes["x"]}')

    # survey one & burst one
    for x in range(chorus_result.sizes["x"]):
        if x % 1000 == 0:
                print(f"{x} samples done for {name}")
       
        variable_survey[x] = np.sum(chorus_result["survey_power"][x,band_range[0]:band_range[1]])
        variable_burst[x] = np.sum(chorus_result["total_power"][x,band_range[0]:band_range[1]])
     
    # save band powers to chorus_result
    chorus_result[f'{name}_power'] = variable_burst
    chorus_result[f'{name}_power'].attrs["description"] = f'integrated power for given (burst, time) coordinate - {name} frequency range'
    
    chorus_result[[f'{name}_survey_power']] = variable_survey
    chorus_result[f'{name}_survey_power'].attrs["description"] = f'survey integrated power for given (burst, time) coordinate - {name} frequency range'



    # save total power to chorus_result
    #power = chorus_result[f'{name}_power']
    # Add attributes (optional)
    #chorus_result["burst_power"].attrs["description"] = "integrated power for given (burst, time) coordinate - full frequency range"

    # find peak power in each band
    for i in range(chorus_result.sizes["x"]):
        #for each frequency band:
        if i % 1000 == 0:
            print(f"{i} peaks done for {name}")
        dt_array = chorus_result[f'{name}_power'].values[i,:]
        dt_array = [np.nan if not np.isfinite(x) else x for x in dt_array]
        dt_array = [np.nan if x==0. else x for x in dt_array]
    
        
        variable_mad[i] = median_abs_deviation(dt_array,nan_policy='omit')
        variable_peak[i] = np.nanmax(dt_array)
        variable_IQR[i] = np.nanmax(dt_array) - np.nanmin(dt_array)
        variable_mean[i] = np.nanmean(dt_array)
        variable_median[i] = np.nanmedian(dt_array)

    
    # save total power to chorus_result
    chorus_result[f'{name}_mad'] = variable_mad
    chorus_result[f'{name}_mad'].attrs["description"] = f'Median averaged diff in int power for given (burst, time) coordinate - {name} frequency range'
    
    chorus_result[[f'{name}_peak']] = variable_peak
    chorus_result[f'{name}_peak'].attrs["description"] = f'Peak int power for given (burst, time) coordinate - {name} frequency range'
    
    chorus_result[f'{name}_IQR'] = variable_IQR
    chorus_result[f'{name}_IQR'].attrs["description"] = f'IQR in int power for given (burst, time) coordinate - {name} frequency range'
    
    chorus_result[[f'{name}_mean']] = variable_mean
    chorus_result[f'{name}_mean'].attrs["description"] = f'Mean in int power for given (burst, time) coordinate - {name} frequency range'
    
    chorus_result[[f'{name}_median']] = variable_mean
    chorus_result[f'{name}_median'].attrs["description"] = f'Median in int power for given (burst, time) coordinate - {name} frequency range'
    
    
    print(f'{name} done!')

chorus_result.to_netcdf('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/all_burst270525.nc')
