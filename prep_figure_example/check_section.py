
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

plt.style.use("ggplot")

def to_dt(nanoseconds_value):

    # Convert nanoseconds to seconds
    seconds_value = nanoseconds_value / 1e9

    # Unix epoch
    epoch = datetime(1970, 1, 1)

    # Convert to datetime
    human_readable_datetime = epoch + timedelta(seconds=seconds_value)

    return human_readable_datetime

""" function that finds indicies of timestamps in list within range [low,high] """
def find_indices_in_range(lst, low, high):
    return [i for i, x in enumerate(lst) if low <= x <= high]

chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24


AE_conds = [(chorus_result["AE"]<100),(chorus_result["AE"]>=100)&(chorus_result["AE"]<300),(chorus_result["AE"]>=300)]
MLAT_cond = (np.abs(chorus_result["MLAT"])>=0.)
#MLT_cond = (np.abs(chorus_result["MLT"])>=6.)

upperband_cond = ((chorus_result["upperband_IQR"]/chorus_result["upperband_median"])>10.)

# set conditions for desired AE range
conditions = (chorus_result[f'lowerband_peak']<1e4)&(chorus_result[f'upperband_peak']<1e4)&(chorus_result[f'lowerband_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'lowerband_survey_power'] !=0) & AE_conds[2] & MLAT_cond & upperband_cond
times =  np.where(conditions,chorus_result['timestamp'][:].values, np.nan)

""" Defining all global variables 
"""
#/data/hpcdata/users/rablack75/power_HPCFlash/out/fig_2015-02-08 01:49:56.185227.png
''' Tell me the date range - start date and end date. Function will also return the number days this range spans '''
#chosen_date = date(year=2016,month=1,day=20)
# set time bounds for burst records
#time_min = datetime(year=2016,month=1,day=20,hour=17,minute=1,second=1)
#time_max = time_min+ timedelta(minutes=60)

for time in times:
    if np.isnan(time):
        continue
    time_dt = to_dt(time)
    print(time_dt)
    chosen_date = time_dt.date()

    """ Create the OMNI dataset """
    omni_dataset = gl.omni_dataset(chosen_date,chosen_date+timedelta(days=2))
    AE, epoch_omni = omni_dataset.omni_stats
    Kp, Dst, epoch_omni_low = omni_dataset.omni_stats_low_res


    # find burst folder
    PSD_folder = '/data/emfisis_burst/wip/rablack75/rablack75/SLURM_out' #'/data/emfisis_burst/wip/rablack75/rablack75/CreateBurst/For_Suman/data/'

    day_files = gl.DataFiles(chosen_date)

    # String version of dates for filename  
    date_string,year,month,day =gl.get_date_string(chosen_date)

    # Get the FFT files
    FFT_files_s = os.path.join(PSD_folder,year,month,'output',day,"PSD_sliding*.cdf")
    FFT_klet  = os.path.join(PSD_folder,year,month,'output',day,"PSD_Kletzing*.cdf")
    FFT_files_s = glob.glob(FFT_files_s)
    FFT_files_k = glob.glob(FFT_klet)

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

    # Getting the LANL data
    lanl_file = h5py.File(day_files.lanl_data)
    lanl_data = gl.AccessLANLAttrs(lanl_file)
    # Getting LANL attributes
    apogee, perigee = lanl_data.apogee_perigee
    Lstar = lanl_data.L_star
    MLT = lanl_data.MLT
    MLAT_N, MLAT_S = lanl_data.MLAT_N_S
    lanl_epoch = lanl_data.epoch

    """ get all densities - including in/out plasmapause flag """
    fpe_uncleaned = density_data.f_pe
    fpe_epoch = density_data.epoch_convert()
    density = density_data.density

    """ the Survey B PSD with trimmings """
    Btotal_reduced = Btotal#fa.BackReduction(Btotal, 'B', False)

    plot_surveyB = {"Bmagnitude": Btotal_reduced, "Frequency": survey_freq, "Epoch": survey_epoch,
                #"fpe": fpe, "fpe_epoch": fpe_epoch,
                "fce": fce, "fce_05": fce_05, "fce_005": fce_005, "f_lhr": f_lhr, "fce_epoch": fce_epoch}

    # save plots to dict
    plots = {"survey":plot_surveyB}
    plots_k ={}
    k = 0    
    for FFT_Sliding,FFT_kletzing in zip(FFT_files_s,FFT_files_k):
        #if k==10:
            #break
        
        FFT_s = pycdf.CDF(FFT_Sliding)
        FFT_k = pycdf.CDF(FFT_kletzing)
        fft_datetime = FFT_s["BurstDatetime"][...]
        if fft_datetime == time_dt:
            #What is the date and time of this burst? Make a folder for that day
            #fft_reduced_s[fft_reduced_s<10**(-9)]=np.nan
            plot_psd= {"PSD": FFT_s["PSD"], "Frequency": FFT_s["Frequencies"], "Time": FFT_s["Time"],
                            "fce": FFT_s["fce"][...],"fce_05": FFT_s["fce_05"][...],  "fce_005": FFT_s["fce_005"][...],
                            "Timestamp":fft_datetime, "flag":'s', "freq lims": [0,1.5*FFT_s["fce"][...]]}
            plot_kpsd= {"PSD": FFT_k["PSD"], "Frequency": FFT_k["Frequencies"], "Time": FFT_k["Time"],
                            "fce": FFT_k["fce"][...],"fce_05": FFT_k["fce_05"][...],  "fce_005": FFT_k["fce_005"][...],
                            "Timestamp":fft_datetime, "flag":'s', "freq lims": [0,1.5*FFT_k["fce"][...]]}
            
            plot_surveyB["Timestamp"] = fft_datetime

            k = k+1
            plots[f"plot{k}" ] = plot_psd
            plots_k[f"plot{k}" ] = plot_kpsd

    

    desired_stamp = plots[f"plot1"]['Timestamp']
    f_lhr_stamp = fp.FindLHR(desired_stamp, fpe_epoch, f_lhr).get_lhr

    """ get power statistics out """

    power_fpath = f'/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/{year}/{month}/power_{date_string}.nc' #f'/data/hpcflash/users/rablack75/power_netCDFs/bug_slurm/{year}/{month}/power_{date_string}.nc'
    #f'/data/emfisis_burst/wip/rablack75/rablack75/CountBurst/CSVs_flashA/bug/{year}/{month}/power_{date_string}.nc'
    load_power = xr.open_dataset(power_fpath)

    """ change the epochs from ns since 1960 to np.datetime64 """

    load_power["time_converted"] = xr.apply_ufunc(
        to_dt,  # The function to apply
        load_power["timestamp"],                 # The xarray DataArray
        vectorize=True            # Whether to vectorize the function (applies to each element)
    )
    
    """ use the desired stamp from psd calculation """

    desired_time = np.datetime64(desired_stamp, 'ns')
    desired_upper = np.datetime64(desired_stamp+timedelta(seconds=6.), 'ns')

    """ turn the conerted times into list, from which we can find the index correspodning to the desired time """
    times = load_power["time_converted"].values
    data_list = [np.datetime64(date) for date in times]
    indices = find_indices_in_range(data_list,desired_time,desired_upper)

    for i in indices:
        if data_list[i] == np.datetime64(desired_stamp):
            #print("woop - located")
            index = i
            #print(index)
            plot_no = 1
            break
        

    """ isolate the power """

    desired_power = load_power["burst_power"].loc[index,:,:]
    load_power["time_converted"].loc[index].values

    # create condition figure
    fig = plt.figure(dpi=300,figsize=(12, 18))
    # GridSpec: 1 rows, 2 columns
    gs = GridSpec(1, 2, figure=fig, width_ratios=[0.9,0.1])
    subfigs = fig.add_subfigure(gs[0])
    axs = subfigs.subplots(4, 1,gridspec_kw={'height_ratios': [2, 2, 2, 2],'hspace': 0.5 })


    name='Accent'
    cmap = mpl.colormaps[name] 
    color_m = cmap.colors
    T_window = 1024*(1/35000)
    duration = 208896*(1/35000)
    n_bins = 406

    t_array = np.linspace(T_window,duration, n_bins).tolist()
    t_array.insert(0,0.)
    t_array = np.array(t_array)
    burst_1 = desired_power

    f05_1 = np.zeros(406)
    f01_05 = np.zeros(406)
    f_all = np.zeros(406)
    exohiss = np.zeros(406)
    for i in range(406):
        f01_05[i] = np.sum(burst_1[1:5,i])
        f05_1[i] = np.sum(burst_1[5:,i])
        f_all[i] = np.sum(burst_1[:,i])
        exohiss[i] = burst_1[0,i]

    non_zero_indices = np.nonzero(f01_05)[0]
    non_zero_indices_upper = np.nonzero(f05_1)[0]
    non_zero_indices_hiss = np.nonzero(exohiss)[0]

    # finding the IQR for each type of wave
    Q1=[]
    Q2=[]
    Q3=[]
    IQR_norm =[]
    for wave in [exohiss,f01_05,f05_1]:

        q1 = np.nanquantile(wave[wave!=0],0.25)
        q2 = np.nanquantile(wave[wave!=0],0.5)
        q3 = np.nanquantile(wave[wave!=0],0.75)

        Q1.append(q1)
        Q2.append(q2)
        Q3.append(q3)

        IQR_norm.append((q3-q1)/q2)
    if IQR_norm[2]<10:
        continue
    axs[3].plot(t_array[non_zero_indices_hiss],exohiss[exohiss != 0],label = f'Exohiss',color=color_m[0])
    axs[3].plot(t_array[non_zero_indices],f01_05[f01_05 != 0],label = f'Lowerband chorus',color=color_m[1])#, linestyle = 'dashed', marker = '*')
    axs[3].plot(t_array[non_zero_indices_upper],f05_1[f05_1 != 0],label = f'Upperband chorus',color=color_m[2])
    #non_zero_indices = np.nonzero(f_all)[0]
    #plt.plot(t_array[non_zero_indices],f_all[f_all != 0],label = f'all')
    #non_zero_indices = np.nonzero(f05_1)[0]
    #plt.plot(t_array[non_zero_indices],f05_1[f05_1 != 0],label = f'upper')
    axs[3].axhline(Q2[0],color = color_m[0],linestyle='dashed')
    axs[3].axhline(Q2[1],color = color_m[1],linestyle='dashed')
    axs[3].axhline(Q2[2],color = color_m[2],linestyle='dashed')

    axs[3].axhline(Q1[0],color = color_m[0],linestyle='dotted')
    axs[3].axhline(Q1[1],color = color_m[1],linestyle='dotted')
    axs[3].axhline(Q2[2],color = color_m[2],linestyle='dotted')

    axs[3].axhline(Q3[0],color = color_m[0],linestyle='dotted')
    axs[3].axhline(Q3[1],color = color_m[1],linestyle='dotted')
    axs[3].axhline(Q3[2],color = color_m[2],linestyle='dotted')

    

    # Add a new axes (invisible) to hold the textbox
    # Coordinates: [left, bottom, width, height] in figure fraction (0 to 1)
    text_ax = subfigs.add_axes([0.95, 0.15, 0.1, 0.1])
    text_ax.axis('off')  # Hide ticks and frame

    # Add text box to the new axes
    text_ax.text(0.1, 0.8, rf'$ \frac{{Q3 - Q1}}{{Q2}} = {(Q3[0] - Q1[0]) / Q2[0]:.2f}$',
                fontsize=12, color=color_m[0])
    text_ax.text(0.1, 0.5, rf'$ \frac{{Q3 - Q1}}{{Q2}} = {(Q3[1] - Q1[1]) / Q2[1]:.2f}$',
                fontsize=12, color=color_m[1])
    text_ax.text(0.1, 0.2, rf'$ \frac{{Q3 - Q1}}{{Q2}} = {(Q3[2] - Q1[2]) / Q2[2]:.2f}$',
                fontsize=12, color=color_m[2])

    axs[3].set_ylabel(r'$Log\ power\ (nT^2)$')
    axs[3].set_xlabel(r'$Time\ (s)$')
    axs[3].set_yscale('log')
    axs[3].legend(loc='upper right')
    axs[3].set_xlim(0,6)
    axs[3].set_title(f'd)) Integrated power for given chorus type', loc='left', fontsize=14)  

    # Large subplot spanning all three columns
    ax1 = axs[0]
    ax1.set_title(f'a) Survey magnetic field spectral density', loc='left', fontsize=14)  
    survey_norm = mcolors.LogNorm(vmin=10**(-9), vmax=10**(-5))

    # plot survey PSD
    plot = plots["survey"]
    survey_img = ax1.pcolormesh(plot["Epoch"], plot["Frequency"], np.transpose(plot["Bmagnitude"]), norm=survey_norm, cmap="viridis", shading='auto')
    # plot gyrofrequencies
    ax1.plot(plot["fce_epoch"],plot["fce"], color = 'white',label = r'$f_{ce}$')
    ax1.plot(plot["fce_epoch"],plot["fce_05"], color = 'white',label = r'$0.5\ f_{ce}$')
    ax1.plot(plot["fce_epoch"],plot["fce_005"], color = 'white',label = r'$0.05\ f_{ce}$')
    # format time x-axis 
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    hours = mdates.HourLocator(interval=1)
    ax1.xaxis.set_minor_locator(hours)
    # set ylims
    ax1.set_ylim(10**1,1.12574*10**4)
    # log scale on y axis
    ax1.set_yscale('log')
    # set labels
    ax1.set_ylabel('Frequency (Hz)',fontsize=8)
    ax1.set_xlabel('UTC',fontsize=8)
    # create legend 
    ax1.legend()
    # Name plot axis
    ax2 = axs[2]
    ax2.set_title(f'c) Burst magnetic field spectral density - 1024-point windows', loc='left', fontsize=14)  
    #cax = axs[1]

    # plot survey PSD
    plot=plots[f"plot1"]
    psd_reduced = plot["PSD"]
    burst_norm = mcolors.LogNorm(vmin=10**(-9), vmax=10**(-5))
    # append final frequency value
    Frequency = plot["Frequency"][:]
    Frequency = np.append(Frequency, plot["Frequency"][-1]+(plot["Frequency"][-1]-plot["Frequency"][-2]))
    burst_samps = ax2.pcolormesh(plot["Time"], np.asarray(Frequency) , np.array(psd_reduced).T, norm=burst_norm, cmap="viridis")

    # plot gyrofrequency lines
    ax2.axhline(plot["fce_05"],0,plot["Time"][-1],color='white',label=r'$0.5\ f_{ce}$')
    ax2.axhline(plot["fce_005"],0,plot["Time"][-1],color='white',label=r'$00.5\ f_{ce}$')
    ax1.axvline(plot["Timestamp"],ymin=0,ymax=np.nanmax(plots["survey"]['fce'])+0.1*np.nanmax(plots["survey"]['fce']),linestyle='dashed',linewidth=2,color='white')

    #for i in range(10):
        #ax2.axhline(freq[i],0,plot["Time"][-1],color='white', linestyle='dashed')

    # add text labels to gyrofreqeuncy lines
    tex1 = r'$0.5\ f_{ce}$'
    ax2.text(4.5, plot["fce_05"]+200, tex1, fontsize=10, va='bottom',color='white')
    tex2 = r'$0.05\ f_{ce}$'
    ax2.text(4.5, plot["fce_005"]+200, tex2, fontsize=10, va='bottom',color='white')

    # set labels
    ax2.set_xlabel('Time (s)',fontsize=8)
    ax2.set_ylabel('Frequency (Hz)',fontsize=8)
    # set ylimits
    ax2.set_ylim(0,2*plot["fce_05"])
    ax2.set_xlim(0, gl.global_constants["Duration"])
    ax2.margins(x=0, y=0)


    # create bounding box for low frequency
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = 0.03       # starting y in data units
    width = 5.95    # width in x-data units
    height = f_lhr_stamp - 0.1*f_lhr_stamp - 0.02 # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[0], facecolor='none', lw=2)
    ax2.add_patch(rect)



    # create bounding box for lower-band
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = f_lhr_stamp+0.2*f_lhr_stamp    # starting y in data units
    width = 5.95    # width in x-data units
    height = plot["fce_05"]-1.3*f_lhr_stamp  # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[1], facecolor='none', lw=2)
    ax2.add_patch(rect)

    # create bounding box for upper-band
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = plot["fce_05"] + 0.2*f_lhr_stamp      # starting y in data units
    width = 5.95    # width in x-data units
    height = plot["fce"]-plot["fce_05"]-0.4*f_lhr_stamp    # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[2], facecolor='none', lw=2)
    ax2.add_patch(rect)

    # create bounding box around where survey product comes from
    ax2.add_patch(plt.Rectangle((0.001, 0.01), 0.078, 0.98, ec='red', fc="none",lw=1,linestyle='dashed',
                            transform=ax2.transAxes))


    # Name plot axis
    ax4 = axs[1]
    ax4.set_title(f'b) Burst magnetic field spectral density - survey 16384-point windows',loc='left', fontsize=14)  
    # plot burst PSD
    plot=plots_k[f"plot1"]
    psd_reduced = plot["PSD"]
    # append final frequency value
    Frequency = plot["Frequency"][:]
    Frequency = np.append(Frequency, plot["Frequency"][-1]+(plot["Frequency"][-1]-plot["Frequency"][-2]))
    burst_samps = ax4.pcolormesh(plot["Time"], np.asarray(Frequency) , np.array(psd_reduced).T, norm=burst_norm, cmap="viridis")

    # set labels
    ax4.set_xlabel('Time (s)',fontsize=8)
    ax4.set_ylabel('Frequency (Hz)',fontsize=8)
    # set ylimits
    ax4.set_ylim(0,2*plot["fce_05"])
    ax4.set_xlim(0, gl.global_constants["Duration"])
    ax4.margins(x=0, y=0)
    # create bounding box for low frequency
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = 0.3*f_lhr_stamp      # starting y in data units
    width = 5.95    # width in x-data units
    height = f_lhr_stamp - 0.4*f_lhr_stamp # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[0], facecolor='none', lw=2)
    ax4.add_patch(rect)



    # create bounding box for lower-band
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = f_lhr_stamp+0.2*f_lhr_stamp    # starting y in data units
    width = 5.95    # width in x-data units
    height = plot["fce_05"]-1.3*f_lhr_stamp  # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[1], facecolor='none', lw=2)
    ax4.add_patch(rect)

    # create bounding box for upper-band
    # Example: rectangle in data coordinates
    x0 = 0.01       # starting x in data units
    y0 = plot["fce_05"] + 0.2*f_lhr_stamp      # starting y in data units
    width = 5.95    # width in x-data units
    height = plot["fce"]-plot["fce_05"]-0.4*f_lhr_stamp    # height in y-data units

    rect = plt.Rectangle((x0, y0), width, height, edgecolor=color_m[2], facecolor='none', lw=2)
    ax4.add_patch(rect)

    # create bounding box around where survey product comes from
    ax4.add_patch(plt.Rectangle((0.001, 0.01), 0.078, 0.98, ec='red', fc="none",lw=1,linestyle='dashed',
                            transform=ax4.transAxes))

    print(ax2.get_ylim(),ax4.get_ylim())

    # do the colorbars:
    # add the colorbar to the leftside
    sfig =subfigs

    # first, for survey power:
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=survey_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.95, 0.75, 0.05, 0.1])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$nT^2$', fontsize=14)
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)

    # Add colorbar for the iqr
    sm = plt.cm.ScalarMappable(cmap='viridis',norm = burst_norm) #norm=iqr_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.95, 0.55, 0.05, 0.1])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$nT^2$', fontsize=14) 
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)

    # Add colorbar for the iqr
    sm = plt.cm.ScalarMappable(cmap='viridis',norm = burst_norm) #norm=iqr_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.95, 0.34, 0.05, 0.1])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01)
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$nT^2$', fontsize=14) 
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)


    plt.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plot_figs_intro/large_upperband_{plot["Timestamp"]}.png')
    print('Saved')