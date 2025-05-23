# Code that sums up power in all frequency bands, saves to total_power variable
# And that finds peak power in all freq bands, and in total power (all frequencies)


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
from scipy.stats import norm,iqr

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




chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/curr_combined/full_13_19.nc')
# Reduce 'chorus_flag' over the 'y' dimension, grouped by 'x'
is_chorus = chorus_result['chorus_flag'].any(dim='y')  # Check for any True in 'y'
chorus_result['isChorus'] = is_chorus

chorus_result = chorus_result.where((chorus_result["isChorus"]),drop=True)

total_power = xr.DataArray(
    data=np.full((chorus_result.sizes["x"],chorus_result.sizes["z"]), np.nan),  # 2D: (x, z)
    dims=("x", "z"),
    name="total_power"
)

print(chorus_result.sizes["x"])
for x in range(chorus_result.sizes["x"]):
    if x % 1000 == 0:
        print(f"{x} done")
    for z in range(chorus_result.sizes["z"]):
        total_power[x,z] = np.sum(chorus_result["burst_power"][x,:,z])

# save total power to chorus_result
chorus_result["total_power"] = total_power
# Add attributes (optional)
chorus_result["total_power"].attrs["description"] = "integrated power for given (burst, time) coordinate - full frequency range"

# find peak power in each sample (all freq)
peak_all_f = np.zeros((chorus_result.sizes["x"]))
for i in range(chorus_result.sizes["x"]):
    # for each frequency band:
    if i % 1000 == 0:
        print(f"{i} done")
    dt_array = chorus_result["total_power"].values[i,:]
    dt_array = [np.nan if not np.isfinite(x) else x for x in dt_array]


    peak_all_f[i] = np.nanmax(dt_array)

total_peak = xr.DataArray(
    data=np.full(chorus_result.sizes["x"], np.nan),  # 2D: (x, z)
    dims="x",
    name="total_power"
)
# save total power to chorus_result
chorus_result["total_peak"] = total_peak
# Add attributes (optional)
chorus_result["total_peak"].attrs["description"] =  "max integrated power for given (burst) coordinate - full frequency range"

# Now find the mean power over the full frequency range

# find mean power in each sample (all freq)

total_mean = xr.DataArray(
    data=np.full(chorus_result.sizes["x"], np.nan),  # 1D: (x)
    dims="x",
    name="total_mean"
)
for i in range(chorus_result.sizes["x"]):
    # for each frequency band:
    if i % 1000 == 0:
        print(f"{i} done")
    dt_array = chorus_result["total_power"].values[i,:]
    dt_array = [np.nan if not np.isfinite(x) else x for x in dt_array]

    total_mean[i] = np.nanmean(dt_array)
    

# save total power to chorus_result
chorus_result["total_mean"] = total_mean
# Add attributes (optional)
chorus_result["total_mean"].attrs["description"] =  "mean integrated power for given (burst) coordinate - full frequency range"



chorus_result.to_netcdf('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/all_burst230525.nc')