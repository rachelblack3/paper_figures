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
from matplotlib.ticker import LogLocator,FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib
matplotlib.use('Agg')

from scipy.stats import ks_2samp
import scipy.stats as stats
# other
import pandas as pd
import string
import roman

plt.style.use("ggplot")


chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus010725_nonlinearthresh.nc')
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

def make_dist(data,bins):

    # Compute the histogram Seaborn is plotting
    counts, bin_edges = np.histogram(
        data,
        bins=bins,
        density=False
    )

    # Convert to probabilities (like stat='probability')
    probabilities = counts / counts.sum()

    # Optional: get bin centers for plotting / analysis
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    return probabilities,bin_centers

def make_subset(conditions, dataset,total_flag,wave,waves,conditions2):

    if total_flag==False:
        survey_vals = np.where(conditions,dataset[f'{name}_survey_power'][:].values, np.nan)
        iqr_vals =  np.where(conditions,dataset[f'{name}_IQR'][:].values, np.nan)
        mad_vals = np.where(conditions,dataset[f'{name}_mad'][:].values, np.nan)
        peak_vals = np.where(conditions,dataset[f'{name}_peak'][:].values, np.nan)
        median_vals = np.where(conditions,dataset[f'{name}_median'][:].values, np.nan)
        ratio_vals = peak_vals/survey_vals
        prob_diff = (peak_vals - survey_vals)/(peak_vals+survey_vals)
        iqr_over_median = iqr_vals/median_vals

        survey_vals =[k for k in survey_vals if k != 0]
        survey_vals = [k for k in survey_vals if k <1e4]
        survey_vals = [k for k in survey_vals if ~np.isnan(k)]

        peak_vals =[k for k in peak_vals if k != 0]
        peak_vals = [k for k in peak_vals if k <1e4]
        peak_vals = [k for k in peak_vals if ~np.isnan(k)]

        ratio_vals =[k for k in ratio_vals if k > 1.]
        #ratio_vals = [k for k in ratio_vals if k <1e4]
        ratio_vals = [k for k in ratio_vals if ~np.isnan(k)]

        iqr_vals =[k for k in iqr_over_median if k != 0]
        #iqr_vals = [k for k in iqr_over_median if k <1e4]
        iqr_vals = [k for k in iqr_vals if ~np.isnan(k)]
    
    elif total_flag==True:
        survey_vals = np.where(conditions[0],dataset[f'{waves[0]}_survey_power'][:].values, np.nan)+ np.where(conditions[1],dataset[f'{waves[1]}_survey_power'][:].values, np.nan)+ np.where(conditions2,dataset[f'{wave}_survey_power'][:].values, np.nan)
        iqr_vals =  np.where(conditions[0],dataset[f'{waves[0]}_IQR'][:].values, np.nan)+np.where(conditions[1],dataset[f'{waves[1]}_IQR'][:].values, np.nan)+ np.where(conditions2,dataset[f'{wave}_IQR'][:].values, np.nan)
        peak_vals = np.where(conditions[0],dataset[f'{waves[0]}_peak'][:].values, np.nan)+np.where(conditions[1],dataset[f'{waves[1]}_peak'][:].values, np.nan)+ np.where(conditions2,dataset[f'{wave}_peak'][:].values, np.nan)
        median_vals = np.where(conditions[0],dataset[f'{waves[0]}_median'][:].values, np.nan)+np.where(conditions[1],dataset[f'{waves[1]}_median'][:].values, np.nan)+ np.where(conditions2,dataset[f'{wave}_median'][:].values, np.nan)
        ratio_vals = peak_vals/survey_vals
        prob_diff = (peak_vals - survey_vals)/(peak_vals+survey_vals)
        iqr_over_median = iqr_vals/median_vals

        survey_vals =[k for k in survey_vals if k != 0]
        survey_vals = [k for k in survey_vals if k <1e4]
        survey_vals = [k for k in survey_vals if ~np.isnan(k)]

        peak_vals =[k for k in peak_vals if k != 0]
        peak_vals = [k for k in peak_vals if k <1e4]
        peak_vals = [k for k in peak_vals if ~np.isnan(k)]

        ratio_vals =[k for k in ratio_vals if k > 1.]
        #ratio_vals = [k for k in ratio_vals if k <1e4]
        ratio_vals = [k for k in ratio_vals if ~np.isnan(k)]

        iqr_vals =[k for k in iqr_over_median if k != 0]
        #iqr_vals = [k for k in iqr_over_median if k <1e4]
        iqr_vals = [k for k in iqr_vals if ~np.isnan(k)]


    return survey_vals,peak_vals,ratio_vals,iqr_vals


def make_total(conditions, dataset):

    survey_vals = np.where(conditions[0],dataset[f'exohiss_band_survey_power'][:].values, np.nan)+ np.where(conditions[1],dataset[f'lowerband_survey_power'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_survey_power'][:].values, np.nan)
    iqr_vals =  np.where(conditions[0],dataset[f'exohiss_band_IQR'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_IQR'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_IQR'][:].values, np.nan)
    peak_vals = np.where(conditions[0],dataset[f'exohiss_band_peak'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_peak'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_peak'][:].values, np.nan)
    median_vals = np.where(conditions[0],dataset[f'exohiss_band_median'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_median'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_median'][:].values, np.nan)
    ratio_vals = peak_vals/survey_vals
    prob_diff = (peak_vals - survey_vals)/(peak_vals+survey_vals)
    iqr_over_median = iqr_vals/median_vals

    
    survey_vals =[k for k in survey_vals if k != 0]
    survey_vals = [k for k in survey_vals if k <1e4]
    survey_vals = [k for k in survey_vals if ~np.isnan(k)]

    peak_vals =[k for k in peak_vals if k != 0]
    peak_vals = [k for k in peak_vals if k <1e4]
    peak_vals = [k for k in peak_vals if ~np.isnan(k)]

    ratio_vals =[k for k in ratio_vals if k > 1.]
    #ratio_vals = [k for k in ratio_vals if k <1e4]
    ratio_vals = [k for k in ratio_vals if ~np.isnan(k)]

    iqr_vals =[k for k in iqr_over_median if k != 0]
    #iqr_vals = [k for k in iqr_over_median if k <1e4]
    iqr_vals = [k for k in iqr_vals if ~np.isnan(k)]

    return survey_vals,peak_vals,ratio_vals,iqr_vals






# define subsets
names = ["exohiss_band", "lowerband", "upperband"]
AE_conds = [(chorus_result["AE"]<100),(chorus_result["AE"]>=100)&(chorus_result["AE"]<300),(chorus_result["AE"]>=300)]
AE_names = [r'$AE<100$',r'$100<=AE< 300$',r'$AE>=300$']
MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)

# make list for D and p_subset
D_subsets = np.zeros((3,3,4))
p_subsets = np.zeros((3,3,4))

# define bins for each statistic
log_start = -10
log_end = 0
log_bins_survey= np.logspace(log_start, log_end, 50)

log_start = -8
log_end = 1
log_bins_peak = np.logspace(log_start, log_end, 50)

log_start = 0
log_end = 4
log_bins_ratio = np.logspace(log_start, log_end, 50)

log_start = -1
log_end = 4
log_bins_iqr = np.logspace(log_start, log_end, 50)

# First the distribution of choice:
for j,name in enumerate(names):
    percent_over_100 = np.zeros(3)
    burst_over_nonlin = np.zeros(3)
    survey_over_nonlin = np.zeros(3)
    iqr_over_10 = np.zeros(3)

    for i, (AE_cond, AE_name) in enumerate(zip(AE_conds, AE_names)):

        print(f'AE condition {i} started...')
        
        # compute your filtered data (independently each time)
        # set conditions for desired AE range
        conditions = (chorus_result[f'{name}_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_survey_power'] !=0) & AE_cond & MLAT_cond

        survey_vals,peak_vals,ratio_vals,iqr_vals = make_subset(conditions,chorus_result,False,False,False,False)
        print(np.nanmin(ratio_vals),np.nanmax(ratio_vals))
        # distribution
        surv_prob, surv_x_values = make_dist(survey_vals,log_bins_survey)
        burst_prob, burst_x_values = make_dist(peak_vals,log_bins_peak)
        ratio_prob, ratio_x_values = make_dist(ratio_vals,log_bins_ratio)
        print(np.nanmin(ratio_x_values),"ratio")
        iqr_prob, iqr_x_values = make_dist(iqr_vals,log_bins_iqr)
        
        # now make the distribution for total-subset
        conditions_exo = (chorus_result[f'exohiss_band_peak']<1e4)&(chorus_result[f'exohiss_band_peak'] !=0)&(chorus_result["Lstar"]>0.) & MLAT_cond
        conditions_lower = (chorus_result[f'lowerband_peak']<1e4)&(chorus_result[f'lowerband_peak'] !=0)&(chorus_result["Lstar"]>0.)& MLAT_cond
        conditions_upper = (chorus_result[f'upperband_peak']<1e4)&(chorus_result[f'upperband_peak'] !=0)&(chorus_result["Lstar"]>0.)& MLAT_cond
        all_conds = [conditions_exo,conditions_lower,conditions_upper]

        conditions2 = (chorus_result[f'{name}_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.) & ~AE_cond & MLAT_cond
        waves = [name for k, name in enumerate(names) if k != j]
        conditions = [conds for k, conds in enumerate(all_conds) if k != j]
        total_survey_vals,total_peak_vals,total_ratio_vals,total_iqr_vals = make_subset(conditions,chorus_result,True,name,waves,conditions2)

        # distribution
        total_surv_prob, total_surv_x_values = make_dist(total_survey_vals,log_bins_survey)
        total_burst_prob, total_burst_x_values = make_dist(total_peak_vals,log_bins_peak)
        total_ratio_prob, total_ratio_x_values = make_dist(total_ratio_vals,log_bins_ratio)
        print(np.nanmin(total_ratio_x_values),"ratio")
        total_iqr_prob, total_iqr_x_values = make_dist(total_iqr_vals,log_bins_iqr)
        
        # full full population 
        full_survey_vals,full_peak_vals,full_ratio_vals,full_iqr_vals = make_total(all_conds,chorus_result)

        # distribution
        full_surv_prob, total_surv_x_values = make_dist(full_survey_vals,log_bins_survey)
        full_burst_prob, total_burst_x_values = make_dist(full_peak_vals,log_bins_peak)
        full_ratio_prob, total_ratio_x_values = make_dist(full_ratio_vals,log_bins_ratio)
        full_iqr_prob, total_iqr_x_values = make_dist(full_iqr_vals,log_bins_iqr)


        def ecdf(data):
            xs = np.sort(data)
            ys = np.arange(1, len(xs)+1)/len(xs)
            return xs, ys

        subset_dists = [survey_vals,peak_vals,ratio_vals,iqr_vals]
        rest_dists = [total_survey_vals,total_peak_vals,total_ratio_vals,total_iqr_vals]
        stat_names=[r'$p_{survey,FFT}$',r'$p_{burst,max}$',r'$p_{survey,FFT}$$p_{burst,max}$',r'$IQR_{norm}$']

        for subset,rest,stat_name in zip(subset_dists,rest_dists,stat_names):
            x_ecdf, Fx = ecdf(np.log10(subset))
            y_ecdf, Fy = ecdf(np.log10(rest))

            # Compute KS statistic
            D, p = ks_2samp(total_iqr_vals, iqr_vals)

            # Plot ECDFs
            plt.figure(figsize=(8,5))
            plt.step(x_ecdf, Fx, where='post', label=f'subset (n={len(subset)})')
            plt.step(y_ecdf, Fy, where='post', label=f'rest-subset (n={len(rest)})')

            # Highlight maximum gap
            idx_max = np.argmax(np.abs(Fx - np.interp(x_ecdf, y_ecdf, Fy)))
            plt.axvline(x_ecdf[idx_max], color='gray', linestyle='--', label='Max gap')
            plt.axhline(Fx[idx_max], color='red', linestyle='--')
            plt.axhline(np.interp(x_ecdf[idx_max], y_ecdf, Fy), color='blue', linestyle='--')
            plt.axhline(0.25, color='green', linestyle='dotted',label='75th percentile')

            if (stat_name == r'$p_{survey,FFT}$'):
                # find where cdf meets nonlinear threshold
                # Find where y crosses the target value
                target = np.log10(0.0225)
                # Find closest index
                y_cross = np.interp(target, x_ecdf, Fx)
                # save to array
                survey_over_nonlin[i]=100*(1-y_cross)
                print(y_cross)
                plt.axvline(np.log10(0.0225), color='red', linestyle='dotted',label='Nonlinear threshold')
                plt.axhline(y_cross, color='red', linestyle='dotted',label='Proportion at nonlinear threshold')
            elif (stat_name == r'$p_{burst,max}$'):
                # find where cdf meets nonlinear threshold
                # Find where y crosses the target value
                target = np.log10(0.0225)
                # Find closest index
                y_cross = np.interp(target, x_ecdf, Fx)
                # save to array
                burst_over_nonlin[i]=100*(1-y_cross)
                print(y_cross)
                plt.axvline(np.log10(0.0225), color='red', linestyle='dotted',label='Nonlinear threshold')
                plt.axhline(y_cross, color='red', linestyle='dotted',label='Proportion at nonlinear threshold')


            if (stat_name == r'$p_{survey,FFT}$$p_{burst,max}$'):
                # find where cdf meets nonlinear threshold
                # Find where y crosses the target value
                target = np.log10(100)
                # Find closest index
                y_cross = np.interp(target, x_ecdf, Fx)

                # save to array
                percent_over_100[i]=100*(1-y_cross)
                print(y_cross,"for the ratio >100")
                plt.axvline(np.log10(0.100), color='red', linestyle='dotted',label='100')
                plt.axhline(y_cross, color='red', linestyle='dotted',label='Proportion at 100')

            if (stat_name == r'$IQR_{norm}$'):
                # find where cdf meets nonlinear threshold
                # Find where y crosses the target value
                target = np.log10(10)
                # Find closest index
                y_cross = np.interp(target, x_ecdf, Fx)

                # save to array
                iqr_over_10[i]=100*(1-y_cross)
                print(y_cross,"for the iqr_norm >10")
                plt.axvline(np.log10(0.100), color='red', linestyle='dotted',label='100')
                plt.axhline(y_cross, color='red', linestyle='dotted',label='Proportion at 100')

            plt.title(f"ECDFs for {stat_name}, {name}, {AE_name}")
            plt.xlabel("Value")
            plt.ylabel("Cumulative Probability")
            plt.legend()
            plt.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/prep_fig4/CDFS/{name}_{AE_name}_{stat_name}.png')

        
    thresholds = [percent_over_100,burst_over_nonlin,survey_over_nonlin,iqr_over_10]

    np.savetxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig4/AE/VanAllenB/LowMLATThresholds_{name}.txt', thresholds, fmt="%.7f") 