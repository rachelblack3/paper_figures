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


# starting with survey, access file that has all survey sampling + chorus occurrence

survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerA_Version20240515_withchorus_pos.csv')

# change all of the timestamps to datetime objects
survey["Timestamp"] = survey["Timestamp"].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
# change MLT to radians for binning in polar coords
survey["MLT_rads"] = 2*np.pi*survey["MLT"]/24

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
samp_mask =  (survey_chorus_smoothed < 0.1*np.nanmedian(survey_chorus_smoothed))
# Print the inf mask
print("\nLow sampling bin mask:")
print(samp_mask)
survey_chorus_smoothed[samp_mask] = np.nan

# CHORUS OCCURRENCE (chorus events/all events)
survey_chorus_occurrence = survey_chorus_smoothed/survey_samp_smoothed



# Now doing burst
burst = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/curr_combined/full_13_18.nc')

# create AE and plasmapause filters
burst_low = burst.where(burst["AE"]<100)
burst_medium= burst.where((burst["AE"]>=100)&(burst["AE"]<300))
data_high =burst.where(burst["AE"]>=300)

# change MLT to radians for binning in polar coords
burst["MLT_rads"] = 2*np.pi*burst["MLT"]/24

# SAMPLING
# make histogram
burst_samp_hist = make_histogram_da(burst)

# remove zeroes
zero_mask = burst_samp_hist == 0
burst_samp_hist[zero_mask]=np.nan

# Kaine's smoothing technique
burst_samp_smoothed = np.zeros((24*10,70))
for i in range(24):
    burst_samp_smoothed[i*10:(i+1)*10,:] = burst_samp_hist[i,:]


# CHORUS EVENTS
# make histogram
burst_chorus_hist = make_histogram_da(burst[(burst["chorus_pos"]==True) & (burst["Plasmapause"]=="Out") & (burst["Lstar"]>2.)])

# remove zeroes
zero_mask = burst_chorus_hist == 0
burst_chorus_hist[zero_mask]=np.nan

# Kaine's smoothing technique
burst_chorus_smoothed = np.zeros((24*10,70))
for i in range(24):
    burst_chorus_smoothed[i*10:(i+1)*10,:] = burst_chorus_hist[i,:]


# CHORUS OCCURRENCE (chorus events/all events)
burst_chorus_occurrence = burst_chorus_smoothed/burst_samp_smoothed


# BURST/SURVEY COMPARISONS
#1. SAMPLING
burst_survey_sampling = burst_samp_smoothed/survey_samp_smoothed

#2. CHORUS OCCURRENCE
burst_survey_chorus = burst_chorus_smoothed/survey_chorus_smoothed

# plot all 6

