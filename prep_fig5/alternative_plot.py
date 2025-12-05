# This code creates the data required to make the bands of power plots (low frequency, lowerband and upperband)
# this time seperated by AE
# CURRENTLY this is using slightly incomplete data as the input whilst the full job is running


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogLocator,FuncFormatter
import matplotlib
matplotlib.use('Agg')
# other
import pandas as pd
import string
import roman

plt.style.use("ggplot")


class LogDivergingNorm(mcolors.Normalize):
    """
    Custom normalization: log-scale diverging colormap centered at 1.
    vmin < 1 < vmax
    """
    def __init__(self, vmin, vmax, vcenter=1):
        super().__init__(vmin, vmax)
        if not (vmin < vcenter < vmax):
            raise ValueError("vmin < vcenter < vmax must hold")
        self.vmin = vmin
        self.vmax = vmax
        self.vcenter = vcenter
        self.log_vmin = np.log10(vmin)
        self.log_vmax = np.log10(vmax)
        self.log_center = np.log10(vcenter)

    def __call__(self, value, clip=None):
        log_val = np.log10(np.clip(value, self.vmin, self.vmax))
        return (log_val - self.log_center) / (self.log_vmax - self.log_vmin) + 0.5


def make_stat_histogram(dataframe,stat):
    
    # define spatial bins: 0.1 L and 1 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)

    power = stat#.tolist()
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
    print(bin_counts,"for ratio")
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


def make_ratio_histogram(dataframe,stat):
    
    # define spatial bins: 0.1 L and 1 MLT bins
    rbins = np.linspace(0,7, 71,endpoint=True)
    abins = np.linspace(0,2*np.pi, 25, endpoint=True)

    power = stat#.tolist()
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
            elif np.isfinite(power[i])==False:
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


def make_dial_plot(ax,histogram,colorbar,title,camp):  
  
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
    pc = ax.pcolormesh(A, R,funct, shading='flat',norm=colorbar_norm,cmap=cmap,edgecolors='face')#vmin=0,vmax=0.1,cmap='jet',edgecolors='face')
    
    
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

def make_log_histogram(vals,log_bins):
    

    nbins = len(log_bins)
    min_bin = log_bins[0]
    max_bin = log_bins[-1]

    hist_array1 = np.zeros(nbins-1)

    for val in vals:
        for i in range(nbins-1):

                    
            if (log_bins[i])<=val<(log_bins[i+1]):
                
                hist_array1[i]+=1
                break

    hist_array1 = hist_array1/np.sum(hist_array1 * np.diff(log_bins))

    return(hist_array1)


# import the chorus dataset
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/all_stats_mean_redidV12_fpefceADDED.nc')
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus270525.nc')
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')

burstA = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/A/full_12_18_withTotalPowerAndSurv.nc')
burstB = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/B/full_12_19_withTotalPowerAndSurv.nc')
chorus_result = xr.concat([burstA, burstB], dim="x")

chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"] 

# Make a 4 x 3 plots (do: survey power, burst peak power, IQR, burst peak/survey)

# Create a figure

#chorus_result = chorus_result.where((chorus_result["isChorus"]),drop=True)

plot_names = list(string.ascii_lowercase)
# Lowercase Roman numerals from 1 to 10
labels = list(string.ascii_lowercase[:18])

fpe_fce_conds = [(chorus_result["fpe_fce"]<5),(chorus_result["fpe_fce"]>=5)&(chorus_result["fpe_fce"]<9),(chorus_result["fpe_fce"]>=9)]
fpe_fce_names = [r'$f_{pe}/f_{ce} < 5$',r'$5 \leq f_{pe}/f_{ce} < 9$',r'$f_{pe}/f_{ce} \geq 9$']
MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)
flags = ['lf','lb','ub']
colors = ['#4477AA','#DDCC77','#CC6677']

linestyles = ['solid','solid','solid']

# define thresholds
ratio_thresh = 24.7
IQR_thresh = 25.4
power_thresh = 2.25e-2
all_thresholds = [power_thresh,power_thresh,ratio_thresh,IQR_thresh]

for (flag,name) in zip(flags,names):
    # define datafile where threshold is being saved
    data_file = f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig4/fpefce/{name}_thresholds'
    ratio_TF = []
    iqr_TF = []
    Burst_TF = []
    Survey_TF = []
    all_thresholds_list = [Survey_TF,Burst_TF,ratio_TF,iqr_TF]

    # create condition figure
    fig = plt.figure(figsize=(18, 16))
    # GridSpec: 1 rows, 4 columns
    gs = GridSpec(1, 6, figure=fig, width_ratios=[0.4,2,2,2,0.6,3])

    # Top row: 4 subfigures
    top_subfigs = [fig.add_subfigure(gs[0, i]) for i in [1,2,3,5]]
    axsRight = top_subfigs[3].subplots(4, 1,) 
    top_subfigs[3].subplots_adjust(hspace=0.5)

    for i,(AE_cond,style,AE_name,col) in enumerate(zip(fpe_fce_conds,linestyles,fpe_fce_names,colors)):

        if i<3:

            axsLeft = top_subfigs[i].subplots(4, 1, subplot_kw={'projection': 'polar'}) 
            top_subfigs[i].subplots_adjust(left=0.05, right=0.85,hspace=0.1)

        print(f'AE condition {i} started...')
        # set conditions for desired AE range
        conditions = (chorus_result[f'{name}_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_surv_power'] !=0) & (chorus_result[f'{flag}_flag']==1)& AE_cond & MLAT_cond

        survey_vals = np.where(conditions,chorus_result[f'{name}_surv_power'][:].values, np.nan)
        iqr_vals =  np.where(conditions,chorus_result[f'{name}_iqr'][:].values, np.nan)
        #mad_vals = np.where(conditions,chorus_result[f'{name}_mad'][:].values, np.nan)
        peak_vals = np.where(conditions,chorus_result[f'{name}_peak'][:].values, np.nan)
        median_vals = np.where(conditions,chorus_result[f'{name}_med'][:].values, np.nan)
        ratio_vals = peak_vals/survey_vals
        prob_diff = (peak_vals - survey_vals)/(peak_vals+survey_vals)
        iqr_over_median = iqr_vals/median_vals
    
        # check the 75th percentiles
        ratio_75 = np.nanpercentile(ratio_vals, 75)
        iqr_75 = np.nanpercentile(iqr_over_median, 75)
        surv_power_75 = np.nanpercentile(survey_vals, 75)
        burst_power_75 = np.nanpercentile(peak_vals, 78)
        all_75 = [surv_power_75,burst_power_75,ratio_75,iqr_75]

        for (value,threshold,list) in zip(all_75,all_thresholds,all_thresholds_list):
            print("The value is:",value," and the threshold is:",threshold)
            if value>threshold:
                list.append(True)
            else:
                list.append(False)
        #ratio_vals = np.where(conditions,ratio, np.nan)
    
        
        # make histograms
        peak_hist,mlt_l_dists = make_stat_histogram(chorus_result,peak_vals)
        peakratio_hist, mlt_l_dists = make_ratio_histogram(chorus_result,ratio_vals)
        surv_hist, mlt_l_dists = make_stat_histogram(chorus_result,survey_vals)
        iqr_hist,mlt_l_dists = make_stat_histogram(chorus_result,iqr_vals)
        #mad_hist,mlt_l_dists = make_stat_histogram(chorus_result,mad_vals)
        iqr_median_hist,mlt_l_dists = make_stat_histogram(chorus_result,iqr_over_median)

        print(f'{name} histograms made...')
        # Kaine's smoothing technique before plotting

        survey_smooth = np.zeros((24*10,70))
        ratio_smooth = np.zeros((24*10,70))
        peak_smooth = np.zeros((24*10,70))
        iqr_smooth = np.zeros((24*10,70))
        mad_smooth= np.zeros((24*10,70))
        iqr_median_smooth = np.zeros((24*10,70))
        for k in range(24):
            survey_smooth[k*10:(k+1)*10,:] = surv_hist[k,:]

            peak_smooth[k*10:(k+1)*10,:] = peak_hist[k,:]

            ratio_smooth[k*10:(k+1)*10,:] = peakratio_hist[k,:]

            iqr_smooth[k*10:(k+1)*10,:] = iqr_hist[k,:]

            #mad_smooth[k*10:(k+1)*10,:] = mad_hist[k,:]

            iqr_median_smooth[k*10:(k+1)*10,:] = iqr_median_hist[k,:]
        

        # define normalisations
        peak_norm = mcolors.LogNorm(vmin=10**(1), vmax=10**(5))
        survey_norm = mcolors.LogNorm(vmin=10**(1), vmax=10**(5))
        ratio_norm = mcolors.LogNorm(vmin=1e1,vmax=1e4)
        #iqr_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
        #iqr_median_norm = mcolors.Normalize(vmin=0., vmax=50.)
        log_frac = iqr_median_smooth
        print(np.nanmin(log_frac))
        print(np.nanmax(log_frac))
        iqr_median_norm = mcolors.LogNorm(vmin=1e-1,vmax=1e2)
        
    
        # save all to list
        norms = [survey_norm,peak_norm,ratio_norm,iqr_median_norm]
        # save all data to list
        datas = [survey_smooth*1e6,peak_smooth*1e6,ratio_smooth,log_frac]
        #norms = [ratio_norm,iqr_median_norm]

        cmaps = ['plasma','plasma','seismic','plasma']

        # now do plot
        print(f'condition {i} plots begun...')
        
        for j,(data,norm,cmap) in enumerate(zip(datas,norms,cmaps)):
            
        
            title = labels[i+4*j] +')'
            axis = axsLeft[j]
            axis.set_facecolor('whitesmoke')
            axis.set_title(title,loc='left',fontsize=14)

            if (j == 0):
                dial_plot = make_dial_plot(axis,data,norm,title,cmap)

                if (i==0):
                    axis.set_xticklabels([str(i) if (i==18 or i == 12) else '' for i in range(0,24)],fontsize=14)
                elif (i==1):
                    axis.set_xticklabels([str(i) if (i == 12) else '' for i in range(0,24)],fontsize=14)
                elif (i==2):
                    axis.set_xticklabels([str(i) if (i==6 or i == 12) else '' for i in range(0,24)],fontsize=14)

            elif (j==1):
                dial_plot = make_dial_plot(axis,data,norm,title,cmap)
                if (i==0):
                    axis.set_xticklabels([str(i) if (i==18) else '' for i in range(0,24)],fontsize=14)
                elif (i==1):
                    axis.set_xticklabels(['' for i in range(0,24)],fontsize=14)
                elif (i==2):
                    axis.set_xticklabels([str(i) if (i==6) else '' for i in range(0,24)],fontsize=14)
            
            elif (j==2):
                dial_plot = make_dial_frac_plot(axis,data,norm,title)
                if (i==0):
                    axis.set_xticklabels([str(i) if (i==18) else '' for i in range(0,24)],fontsize=14)
                elif (i==1):
                    axis.set_xticklabels(['' for i in range(0,24)],fontsize=14)
                elif (i==2):
                    axis.set_xticklabels([str(i) if (i==6) else '' for i in range(0,24)],fontsize=14)
            
            elif (j == 3):
            
                dial_plot = make_dial_plot(axis,data,norm,title,cmap)
                if (i==0):
                    axis.set_xticklabels([str(i) if (i==18 or i == 0) else '' for i in range(0,24)],fontsize=14)
                elif (i==1):
                    axis.set_xticklabels([str(i) if (i == 0) else '' for i in range(0,24)],fontsize=14)
                elif (i==2):
                    axis.set_xticklabels([str(i) if (i==6 or i == 0) else '' for i in range(0,24)],fontsize=14)
            
            
                

            #box = axis.get_position()
            #axis.set_position([box.x0 - 0.05, box.y0, box.width, box.height])

        
        survey_vals =[k for k in survey_vals if k != 0]
        survey_vals = [k for k in survey_vals if k <1e4]
        survey_vals = [k*1e6 for k in survey_vals]
        log_start =-2 #-10
        log_end = 6
        log_bins = np.logspace(log_start, log_end, 50)
        #log_bins = np.logspace(np.log10(np.nanmin(survey_vals)), np.log10(np.nanmax(survey_vals)), 50)
        print(log_bins)
        hist_survey = make_log_histogram(survey_vals,log_bins)
        hist_survey_bins = log_bins
        sns.histplot(data=survey_vals,bins=log_bins,color=col,label=AE_name, kde=False,stat='probability',fill=False,element="step",ax=axsRight[0],linestyle=style)

        peak_vals =[k for k in peak_vals if k != 0]
        peak_vals = [k for k in peak_vals if k <1e4]
        peak_vals = [k*1e6 for k in peak_vals]
        log_start = -1
        log_end = 7
        log_bins = np.logspace(log_start, log_end, 50)
        hist_peak = make_log_histogram(peak_vals,log_bins)
        hist_peak_bins = log_bins
        sns.histplot(data=peak_vals,bins=log_bins,color=col,label=AE_name, kde=False,stat='probability',fill=False,element="step",ax=axsRight[1],linestyle=style)

        
        ratio_vals = [i for i in ratio_vals if i >10**(0)]
        ratio_vals = [i for i in ratio_vals if i <1e10]
        #log_bins = np.logspace(np.log10(np.nanmin(ratio_vals)), np.log10(np.nanmax(ratio_vals)), 50)
        log_start = 0
        log_end = 4
        log_bins = np.logspace(log_start, log_end, 50)
        print(log_bins)
        sns.histplot(data=ratio_vals,bins=log_bins,color=col,label=AE_name, kde=False,stat='probability',element="step",fill=False,ax=axsRight[2],linestyle=style)

        log_start = -1
        log_end = 3
        log_bins = np.logspace(log_start, log_end, 50)
        hist_iqr = make_log_histogram(iqr_over_median,log_bins)
        hist_iqr_bins = log_bins
        sns.histplot(data=iqr_over_median,bins=log_bins,color=col,label=AE_name, kde=False,stat='probability',element="step",fill=False,ax=axsRight[3],linestyle=style)

    # save datafile
    df = pd.DataFrame({
    "SurvPow": Survey_TF,
    "BurstPow": Burst_TF,
    "Ratio": ratio_TF,
    "IQR": iqr_TF
    })

    df.to_csv(data_file, index=False)

    title = 'd)'
    axsRight[0].set_title(title,loc='left',fontsize=14)
    axsRight[0].set_xscale('log')
    axsRight[0].axvline(0.0225*10**6,color='black',linestyle='dashed',label='Nonlinear threshold')
    sns.set(font_scale=0.8)
    axsRight[0].set_ylabel(r'Probabillity',fontsize=14)
    axsRight[0].set_xlabel(r'$P_{survey, FFT}\ (pT^2)$',fontsize=14)
    axsRight[0].legend(loc='upper right', bbox_to_anchor=(1, 1))
    axsRight[0].set_facecolor('whitesmoke')
    axsRight[0].set_xlim(10**(-3),10**(8))
    axsRight[0].xaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0], numticks=20))
    axsRight[0].xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=100))

        # Custom formatter: only label every other decade
    def skip_every_other(x, pos):
        exponent = int(np.log10(x))
        if exponent % 2 == 0:  # label only even decades
            return r"$10^{{{}}}$".format(exponent)
        else:
            return ""  # hide odd decades

    axsRight[0].xaxis.set_major_formatter(FuncFormatter(skip_every_other))


    title = 'h)'
    axsRight[1].set_title(title,loc='left',fontsize=14)
    axsRight[1].set_xscale('log')
    sns.set(font_scale=0.8)
    axsRight[1].set_ylabel(r'Probabillity',fontsize=14)
    axsRight[1].set_xlabel(r'$P_{burst, max}\ (pT^2)$',fontsize=14)
    axsRight[1].axvline(0.0225*10**6,color='black',linestyle='dashed',label='Nonlinear threshold')
    #axsRight[1].legend(loc="upper right")
    axsRight[1].set_facecolor('whitesmoke')
    axsRight[1].set_xlim(10**(-3),10**(8))
    axsRight[1].xaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0], numticks=20))
    axsRight[1].xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=100))
    axsRight[1].xaxis.set_major_formatter(FuncFormatter(skip_every_other))


    title = 'l)'
    axsRight[2].set_title(title,loc='left',fontsize=14)
    axsRight[2].set_xscale('log')
    sns.set(font_scale=0.8)
    axsRight[2].set_ylabel(r'Probabillity',fontsize=14)
    axsRight[2].set_xlabel(r'$P_{burst, max}$/$P_{survey, FFT}$',fontsize=14)
    #axsRight[2].legend(loc="upper right")
    axsRight[2].set_facecolor('whitesmoke')
    axsRight[2].set_xlim(10**(0),10**(4))
    axsRight[2].tick_params(axis='both', which='both', direction='out', length=2)

    title = 'p)'
    axsRight[3].set_title(title,loc='left',fontsize=14)
    axsRight[3].set_xscale('log')
    sns.set(font_scale=0.8)
    axsRight[3].set_ylabel(r'Probabillity',fontsize=14)
    axsRight[3].set_xlabel(r'$IQR_{norm}$',fontsize=14)
   # axsRight[3].legend(loc="upper right")
    axsRight[3].set_facecolor('whitesmoke')
    axsRight[3].set_xlim(10**(-1),10**(3))
    axsRight[3].tick_params(axis='both', which='both', direction='out', length=2)



    top_subfigs[0].suptitle(r'$f_{pe}/f_{ce} <5$', x=0.0, ha='left',fontsize=20)
    top_subfigs[1].suptitle(r'$5\leq f_{pe}/f_{ce} <9$',x=0.0, ha='left', fontsize=20)
    top_subfigs[2].suptitle(r'$f_{pe}/f_{ce} \geq 9$', x=0.0, ha='left',fontsize=20)
    top_subfigs[3].suptitle(r'$probabillity\ distribtuion\ of\ all\ cases$', x=0.0, ha='left',fontsize=20)

    # add the colorbar to the leftside
    sfig =top_subfigs[2]
    # first, for power:
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=survey_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.98, 0.6, 0.05, 0.2])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01,ticks=[1e1, 1e3, 1e5])
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$pT^2$', fontsize=14)
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)

    # add the colorbar to the leftside
    sfig =top_subfigs[2]
    # first, for power:
    sm = plt.cm.ScalarMappable(cmap='Greens', norm=ratio_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.98, 0.3, 0.05, 0.2])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01,ticks=[10, 100, 1000])
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    # Remove grey background and optional frame
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)

    # Add colorbar for the iqr
    sm = plt.cm.ScalarMappable(cmap='plasma',norm = iqr_median_norm) #norm=iqr_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.98, 0.1, 0.05, 0.2])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.01,ticks=[10**(0), 10**(2)])
    cbar.set_ticklabels([r'$10^0$',r'$10^2$'])
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar_ax.set_facecolor('none')            # makes it transparent
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)

    # Add side labels for each row
    fig.text(0.01, 0.8, r'$P_{survey, FFT}$', va='center', ha='center', fontsize=20, rotation='vertical')
    fig.text(0.01, 0.6, r'$P_{burst,max}$', va='center', ha='center', fontsize=20, rotation='vertical')
    fig.text(0.01, 0.4,r'$P_{burst,max}$/$P_{survey, FFT}$' , va='center', ha='center', fontsize=20, rotation='vertical')
    fig.text(0.01, 0.2, r'$IQR_{norm}$', va='center', ha='center', fontsize=20, rotation='vertical')

    print(f'{name} done')
            
    fig.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig5/combinedAB/low/{name}_31225_pT2.png')