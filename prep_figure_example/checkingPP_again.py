
# numerical essentials 
import numpy as np

# for plotting
from matplotlib import pyplot as plt, animation
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from metpy.plots import add_timestamp
import os
import pandas as pd
import xarray as xr

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable

# for cdf reading
from spacepy import toolbox
#spacepy.toolbox.update(leapsecs=True)
from spacepy import pycdf
import h5py

# other numerical tools

import os
import glob
from datetime import datetime,timedelta,date

import global_use as gl
import find_properties as fp



""" Defining all global variables 
"""
#/data/hpcdata/users/rablack75/power_HPCFlash/out/fig_2015-02-08 01:49:56.185227.png
''' Tell me the date range - start date and end date. Function will also return the number days this range spans '''
chosen_date = date(year=2016,month=1,day=21)

""" Create the OMNI dataset """
omni_dataset = gl.omni_dataset(chosen_date,chosen_date+timedelta(days=2))
AE, epoch_omni = omni_dataset.omni_stats
Kp, Dst, epoch_omni_low = omni_dataset.omni_stats_low_res

day_files = gl.DataFiles(chosen_date)

# String version of dates for filename  
date_string,year,month,day =gl.get_date_string(chosen_date)


""" Getting survey file and accessing survey frequencies, epoch and magnitude """
survey_file = pycdf.CDF(day_files.survey_data)
survey_data = gl.AccessSurveyAttrs(survey_file)

""" get properties """
survey_freq = survey_data.frequency
survey_epoch = survey_data.epoch_convert()
survey_bins = survey_data.bin_edges
Btotal = survey_data.Bmagnitude
Etotal = survey_data.Emagnitude
#E_reduced = fa.BackReduction(Etotal, 'E', False)

""" density/magnetometer/LANL """
# get the density, LANL, and magentometer data
# Getting the magnetometer data
mag_file = pycdf.CDF(day_files.magnetic_data)
mag_data = gl.AccessL3Attrs(mag_file)
fce, fce_05, fce_005,f_lhr = mag_data.f_ce
fce_epoch = mag_data.epoch_convert()

if day_files.l4_data()[1] == False:
    print("woah")
elif day_files.l4_data()[1] == True:
    print("yaoh")

# Getting the density data
density_file = pycdf.CDF(day_files.l4_data()[0])
density_data = gl.AccessL4Attrs(density_file)

""" Getting the LANL data """
lanl_file = h5py.File(day_files.lanl_data)
lanl_data = gl.AccessLANLAttrs(lanl_file)
# Getting LANL attributes
apogee, perigee = lanl_data.apogee_perigee
Lstar = lanl_data.L_star
MLT = lanl_data.MLT
MLAT_N, MLAT_S = lanl_data.MLAT_N_S
lanl_epoch = lanl_data.epoch

# Find what the first and last half orbits begin with 
all_times = np.concatenate((apogee,perigee))

    # Creating corresponding labels
labels = np.array(['a'] * len(apogee) + ['p'] * len(perigee))

# Sort dates and labels together
sorted_pairs = sorted(zip(all_times, labels))

# Unzip the sorted pairs
all_times, labels = zip(*sorted_pairs)

# Finding the label for the maximum value
max_index,min_index= np.argmax(all_times),np.argmin(all_times)  # Index of the max value
max_label,min_label= labels[max_index],labels[min_index]  # Corresponding label
print(max_index)


fpe_uncleaned = density_data.f_pe
fpe_epoch = density_data.epoch_convert()
density_uncleaned= density_data.density
# clean fpe array
fpe = np.zeros((len(density_uncleaned)))
fpe[:] = np.nan
density = np.zeros((len(density_uncleaned)))
density[:] = np.nan


# Set all of the rogue density/fpe values as np.nan
for i in range(len(density)):

    if fpe_uncleaned[i] <0.:
        density[i] = np.nan
        fpe[i] = np.nan

    else:
        density[i] = density_uncleaned[i]
        fpe[i] = fpe_uncleaned[i]

density_of_interest=[]
time_of_interest=[]
lstar_of_interest=[]
for gee_index in range(len(all_times)-1):
    if labels[gee_index] == 'a':
        for j,time in enumerate(fpe_epoch):
            if all_times[2]<time<all_times[3]:
                print(all_times[2])
                print("On inward journey")
                findLstar= fp.FindLANLFeatures(fpe_epoch[j], lanl_epoch, MLT, MLAT_N, MLAT_S, Lstar)
                Lstar_stamp = findLstar.get_Lstar
                lstar_of_interest.append(10*(6.6/Lstar_stamp)**4)
                density_of_interest.append(density[j])
                print(j,density[j])
        break


plt.plot(density_of_interest)
plt.plot(lstar_of_interest)
plt.axhline(50)
plt.yscale('log')
plt.savefig("density_inbound.png")