import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import xarray as xr
import pandas as pd
# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches

from datetime import datetime

# functions for making the power dial plots

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


# also read in the data from the last figure on sampling: low MLAT and high MLAT
#read_data = np.loadtxt('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/LowMLATdataV2.txt')
#read_data = read_data.reshape((6,24*10,70))
names = ["Survey_sampling","Burst_sampling","Burst_survey_sampling","Survey_chorus","Burst_chorus","Burst_survey_chorus"]
arrays = [np.loadtxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig1/VanAllenB/LowMLATdataV1_{name}.txt') for name in names]
read_data = np.stack(arrays)  # shape (6,240,70)

# read survey and burst chorus occurrence
survey_chorus_occurrence = read_data[3,:,:]/read_data[0,:,:]
burst_chorus_occurrence = read_data[4,:,:]/read_data[1,:,:]

print("reading survey...")
survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerB_Version20250827.csv')
# do the MLT fix
survey["MLT_rads"] = 2*np.pi*survey["MLT"]/24
# change all of the timestamps to datetime objects
survey["Timestamp"] = survey["Timestamp"].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))

# keep only between 2013-2019 for now
time_condition = (survey["Timestamp"] > datetime(year=2012,month=12,day=31)) & (survey["Timestamp"] < datetime(year=2019,month=12,day=30))
survey = survey[time_condition]
# find the chorus in survey 
survey_chorus = survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>2.) & (np.abs(survey["MLAT"])<6.)]

# get the power for these survey events and put into a hist
survey_power_hist,mlt_l_dists = make_stat_histogram(survey_chorus,survey_chorus["Survey_integral"].values)
print("are these lens the same?",len(survey_chorus["Survey_integral"]),len(survey["Survey_integral"]))

survey_power_smooth = np.zeros((24*10,70))
for i in range(24):
    survey_power_smooth[i*10:(i+1)*10,:] = survey_power_hist[i,:]
    
# also remove bins for which there was low sampling for survey, which would skew the statistics. We have chosen 1000 samples as lowest
survey_samp_mask = np.isnan(survey_chorus_occurrence)
survey_power_smooth[survey_samp_mask] = np.nan

# Multiply the power by the chorus occurrence, so that we can get the effective power (?maybe different name)
relative_survey_power = survey_power_smooth*survey_chorus_occurrence


# Now onto burst power...
#burst_chor = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/all_stats_mean_redidV12_fpefceADDED.nc')
#burst_chor = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/curr_combined/full_13_19.nc')
burst_chor = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/all_burstB270725.nc')
is_chorus = burst_chor['chorus_flag'].any(dim='y')  # Check for any True in 'y'
burst_chor['isChorus'] = is_chorus

chorus_result = burst_chor.where((burst_chor["isChorus"]==True))

chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

mean_burst_int = xr.DataArray(
    data=np.full((chorus_result.sizes["x"]), np.nan),  # 2D: (x, z)
    dims=("x"),
    name="mean_burst_int"
)

for x in range(chorus_result.sizes["x"]):
    if x % 1000 == 0:
        print(f"{x} done")
    mean_burst_int[x] = np.nanmean(chorus_result["total_power"][x,:])

# save total power to chorus_result
chorus_result["mean_burst_int"] = mean_burst_int
chorus_result.to_netcdf('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/burstB270725_justChorus.nc')
# define latitude condition
conditions = np.abs(chorus_result["MLAT"])<6.

# seperate out 'mean power' in each 6 seconds!
mean_vals = np.where(conditions,chorus_result["mean_burst_int"][:].values, np.nan)

# get corresponding survey vals: # survey one & burst one
variable_survey = xr.DataArray(
    data=np.full(chorus_result.sizes["x"], np.nan),  # 2D: (x, z)
    dims=("x"),
    name="survey_power"
)
for x in range(chorus_result.sizes["x"]):
    if x % 1000 == 0:
            print(f"{x} samples done")
    
    variable_survey[x] = np.sum(chorus_result["survey_power"][x,:])

chorus_result["tot_survey"] = variable_survey

survey_vals = np.where(conditions,chorus_result["tot_survey"][:].values, np.nan)

# make histogram
mean_chorus_hist,mlt_l_dists = make_stat_histogram(chorus_result,mean_vals)
survey_for_burst_chorus_hist,mlt_l_dists = make_stat_histogram(chorus_result,survey_vals)

survey_for_burst_smooth = np.zeros((24*10,70))
burst_power_smooth = np.zeros((24*10,70))
for i in range(24):
    burst_power_smooth[i*10:(i+1)*10,:] = mean_chorus_hist[i,:]
    survey_for_burst_smooth[i*10:(i+1)*10,:] = survey_for_burst_chorus_hist[i,:]
    
# also remove bins for which there was low sampling for survey, which would skew the statistics. We have chosen 1000 samples as lowest
burst_samp_mask = np.isnan(burst_chorus_occurrence)
burst_power_smooth[burst_samp_mask] = np.nan
survey_for_burst_smooth[burst_samp_mask] = np.nan

# Multiply the power by the chorus occurrence, so that we can get the effective power (?maybe different name)
relative_burst_power = burst_power_smooth*burst_chorus_occurrence

dists = [survey_power_smooth,relative_survey_power,burst_power_smooth,relative_burst_power,survey_for_burst_smooth]
names = ["Survey_power_just_chorus_events","Survey_power_all_events","burst_power_just_chorus_events","burst_power_all_events","survey_for_burst"]

# save all data to a file so we can fix the plotting faster 
# Write the array to disk

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for dist,name in zip(dists,names):


# The formatting string indicates that I'm writing out
# the values in left-justified columns 7 characters in width
# with 2 decimal places. 

    np.savetxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig2/VanAllenB/LowMLATdataV1_{name}.txt', dist) 

    # Writing out a break to indicate different dists...
    #outfile.write(f'{name}')