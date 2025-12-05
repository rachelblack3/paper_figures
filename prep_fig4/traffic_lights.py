# define subsets
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
from scipy.stats import ks_2samp
import scipy.stats as stats
# other
import pandas as pd
import string
import roman

plt.style.use("ggplot")

def make_subset(conditions, dataset,total_flag,wave,waves,conditions2,all_conds):

    if total_flag==False:
        survey_vals = np.where(conditions,dataset[f'{name}_surv_power'][:].values, np.nan)
        iqr_vals =  np.where(conditions,dataset[f'{name}_iqr'][:].values, np.nan)
        #mad_vals = np.where(conditions,dataset[f'{name}_mad'][:].values, np.nan)
        peak_vals = np.where(conditions,dataset[f'{name}_peak'][:].values, np.nan)
        median_vals = np.where(conditions,dataset[f'{name}_med'][:].values, np.nan)
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
        
        survey_vals = np.where(all_conds[0],dataset[f'{waves[0]}_surv_power'][:].values, np.nan)+ np.where(all_conds[1],dataset[f'{waves[1]}_surv_power'][:].values, np.nan)+ np.where(~conditions,dataset[f'{wave}_surv_power'][:].values, np.nan)
        iqr_vals =  np.where(all_conds[0],dataset[f'{waves[0]}_iqr'][:].values, np.nan)+np.where(all_conds[1],dataset[f'{waves[1]}_iqr'][:].values, np.nan)+ np.where(~conditions,dataset[f'{wave}_iqr'][:].values, np.nan)
        peak_vals = np.where(all_conds[0],dataset[f'{waves[0]}_peak'][:].values, np.nan)+np.where(all_conds[1],dataset[f'{waves[1]}_peak'][:].values, np.nan)+ np.where(~conditions,dataset[f'{wave}_peak'][:].values, np.nan)
        median_vals = np.where(all_conds[0],dataset[f'{waves[0]}_med'][:].values, np.nan)+np.where(all_conds[1],dataset[f'{waves[1]}_med'][:].values, np.nan)+ np.where(~conditions,dataset[f'{wave}_med'][:].values, np.nan)
        
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

    survey_vals = np.where(conditions[0],dataset[f'exohiss_band_surv_power'][:].values, np.nan)+ np.where(conditions[1],dataset[f'lowerband_surv_power'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_surv_power'][:].values, np.nan)
    iqr_vals =  np.where(conditions[0],dataset[f'exohiss_band_iqr'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_iqr'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_iqr'][:].values, np.nan)
    peak_vals = np.where(conditions[0],dataset[f'exohiss_band_peak'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_peak'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_peak'][:].values, np.nan)
    median_vals = np.where(conditions[0],dataset[f'exohiss_band_med'][:].values, np.nan)+np.where(conditions[1],dataset[f'lowerband_med'][:].values, np.nan)+ np.where(conditions[2],dataset[f'upperband_med'][:].values, np.nan)
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

# 2️⃣ Function to compute bootstrap CIs
def bootstrap_ci(data1, data2, quantile=0.5, n_boot=5000, ci=95):
    """Bootstrap CI for difference in a quantile (default median = 0.5)."""
    boot_diffs = []
    n1, n2 = len(data1), len(data2)

    for _ in range(n_boot):
        resample1 = np.random.choice(data1, n1, replace=True)
        resample2 = np.random.choice(data2, n2, replace=True)
        q1 = np.quantile(resample1, quantile)
        q2 = np.quantile(resample2, quantile)
        boot_diffs.append(q1 - q2)

    boot_diffs = np.array(boot_diffs)
    lower = np.percentile(boot_diffs, (100 - ci) / 2)
    upper = np.percentile(boot_diffs, 100 - (100 - ci) / 2)
    return boot_diffs, (lower, upper)


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

burstA = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/A/full_12_18_withTotalPowerAndSurv.nc')
burstB = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/B/full_12_19_withTotalPowerAndSurv.nc')
chorus_result = xr.concat([burstA, burstB], dim="x")

chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24



names = ["exohiss_band", "lowerband", "upperband"]
AE_conds = [(chorus_result["AE"]<100),(chorus_result["AE"]>=100)&(chorus_result["AE"]<300),(chorus_result["AE"]>=300)]
AE_names = [r'$AE<100$',r'$100<=AE<300$',r'$AE>=300$']
MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)

# make list for D and p_subset
D_subsets = np.zeros((3,3,4))
p_subsets = np.zeros((3,3,4))

log_start = -7.5
log_end = 0
log_bins_survey = np.logspace(log_start, log_end, 50)


log_start = -5.5
log_end = 1
log_bins_peak = np.logspace(log_start, log_end, 50)

log_start = 0
log_end = 4
log_bins_ratio = np.logspace(log_start, log_end, 50)

log_start = -1
log_end = 2
log_bins_iqr = np.logspace(log_start, log_end, 50)

flags=['lf','lb','ub']

lower_bound = np.zeros((3,3,4,3)) # shpae is freq band, ae range, stat, percentile
upper_bound = np.zeros((3,3,4,3)) # shpae is freq band, ae range, stat, percentile


# First the distribution of choice:
for j,(flag,name) in enumerate(zip(flags,names)):
    print(f'freq band {name} started...')
    for i, (AE_cond, AE_name) in enumerate(zip(AE_conds, AE_names)):

        print(f'AE condition {i} started...')
        
        # compute your filtered data (independently each time)
        # set conditions for desired AE range
        conditions = (chorus_result[f'{flag}_flag'] == 1) & MLAT_cond & AE_cond

        survey_vals,peak_vals,ratio_vals,iqr_vals = make_subset(conditions,chorus_result,False,False,False,False,False)
        print(np.nanmin(ratio_vals),np.nanmax(ratio_vals))
        # distribution
        surv_prob, surv_x_values = make_dist(survey_vals,log_bins_survey)
        burst_prob, burst_x_values = make_dist(peak_vals,log_bins_peak)
        ratio_prob, ratio_x_values = make_dist(ratio_vals,log_bins_ratio)
        print(np.nanmin(ratio_x_values),"ratio")
        iqr_prob, iqr_x_values = make_dist(iqr_vals,log_bins_iqr)
        
        # now make the distribution for total-subset
        
        conditions_exo = (chorus_result['lf_flag'] == 1) & MLAT_cond
        conditions_lower = (chorus_result['lb_flag'] == 1)& MLAT_cond
        conditions_upper = (chorus_result['ub_flag'] == 1)& MLAT_cond
        all_conds = [conditions_exo,conditions_lower,conditions_upper]

        conditions2 =  ~(chorus_result[f'{flag}_flag']==1 & AE_cond & MLAT_cond)         # except when x & parameter are True

        waves = [name for k, name in enumerate(names) if k != j]
        conditions = ((chorus_result[f'{flag}_flag'] == 1) & MLAT_cond & AE_cond)

        
        total_survey_vals,total_peak_vals,total_ratio_vals,total_iqr_vals = make_subset(conditions,chorus_result,True,name,waves,conditions2,all_conds)

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

        # just save subset burst powers

        # 

        for k, (subset,rest,stat_name) in enumerate(zip(subset_dists,rest_dists,stat_names)):

            
            # 3️⃣ Compute bootstrap CIs for median (50th percentile)
            l = 0 # the percentile index in lower_bound and upper_bound array
            boot_25, ci_25 = bootstrap_ci(np.log10(subset), np.log10(rest), quantile=0.25)
            print(f"25th percentile difference CI (95%): {ci_25}")
            lower_bound[j,i,k,l] = ci_25[0]
            upper_bound[j,i,k,l] = ci_25[1]

            l = 1
            # 3️⃣ Compute bootstrap CIs for median (50th percentile)
            boot_median, ci_median = bootstrap_ci(np.log10(subset), np.log10(rest), quantile=0.5)
            print(f"Median difference CI (95%): {ci_median}")
            lower_bound[j,i,k,l] = ci_median[0]
            upper_bound[j,i,k,l] = ci_median[1]

            l = 2
            # 4️⃣ Compute bootstrap CIs for 75th percentile
            boot_p75, ci_p75 = bootstrap_ci(np.log10(subset), np.log10(rest), quantile=0.75)
            print(f"75th percentile difference CI (95%): {ci_p75}")
            lower_bound[j,i,k,l] = ci_p75[0]
            upper_bound[j,i,k,l] = ci_p75[1]

            l = 3
            # 4️⃣ Compute bootstrap CIs for 75th percentile
            boot_p90, ci_p90 = bootstrap_ci(np.log10(subset), np.log10(rest), quantile=0.9)
            print(f"75th percentile difference CI (95%): {ci_p90}")
            lower_bound[j,i,k,l] = ci_p90[0]
            upper_bound[j,i,k,l] = ci_p90[1]


            # 5️⃣ Plot bootstrap distributions
            plt.figure(figsize=(8,4))
            plt.hist(boot_median, bins=40, alpha=0.6, label='Median diffs')
            plt.hist(boot_p75, bins=40, alpha=0.6, label='75th percentile diffs')
            plt.axvline(0, color='k', linestyle='--')
            plt.legend()
            plt.title(f'Bootstrap distributions of quantile differences (subset - the rest) for {stat_name} and {AE_name}')
            plt.xlabel("Difference")
            plt.ylabel("Frequency")
            plt.show()





# shpae is freq band, ae range, stat, percentile
print("For the burst power")
print(lower_bound[0,:,1,2])
print(lower_bound[1,:,1,2])
print(lower_bound[2,:,1,2])

print("For the survey power")
print(lower_bound[0,:,0,2])
print(lower_bound[1,:,0,2])
print(lower_bound[2,:,0,2])

# magnitude of differene
threshold = 3
traffic_AE = np.zeros((3,3))
for j in range(3):
    for k in range(3):
        if 10**upper_bound[j,k,2,2]>threshold:
            traffic_AE[j,k]+=1
        if 10**upper_bound[j,k,3,2]>threshold:
            traffic_AE[j,k]+=1

# add nonlinear thresh:
traffic_AE[1,2]+=1
traffic_AE[1,1]+=1


# make traffic light plots
names = ["exohiss_band", "lowerband", "upperband"]
conditions = ["AE","fpefce"]
labels = [[r'$AE < 100$',r'$100 \leq AE < 300$',r'$AE \geq 300$'],[r'$f_{pe}/f_{ce} < 5$',r'$5 \leq f_{pe}/f_{ce} < 9$',r'$f_{pe}/f_{ce} \geq 9$']]
from matplotlib import colors


# Plot
plt.style.use("ggplot")
fig,axes = plt.subplots(1,2,figsize=(10,5))
ax = axes

array = traffic_AE
   
color_table = np.full((3,3),'black',dtype="<U10")

for j,name in enumerate(names):

    for k in range(3):

        positives = array[j, k]
        print(positives)

        if positives==0:
            color_table[j,k] = 'g'
            print(j,k,'green')
        if positives==1:
            color_table[j,k] = 'yellow'
            print(j,k,'yellow')
        if positives==2:
            color_table[j,k] = 'orange'
            print(j,k,'orange')
        if positives==3:
            color_table[j,k] = 'r'
            print(j,k,'red')
        
# Convert color names to RGBA
rgba_array = np.vectorize(colors.to_rgba)(color_table)

# Reshape into (3, 3, 4) for imshow
rgba_array = np.moveaxis(rgba_array, 0, -1)


ax.imshow(rgba_array, aspect="equal",interpolation='none', extent=(0, 3, 0, 3))

# Set tick labels at cell centers
ax.set_xticks([])
ax.set_yticks([])


#ax.set_xticklabels(labels[i])
#axes[0].set_yticklabels(["Low frequency","Lower-band","Upper-band"])
# Remove the actual tick marks (no lines)

chorus = ["Upper-band","Lower-band","Low frequency"]
# Add text labels centered in cells
for m in range(3):
    ax.text(m+0.5, -0.3, labels[i][m],
                ha='center', va='center', color='black', fontsize=12
                )

for j in range(3):
    axes[0].text(-0.6, j+0.5, chorus[j],
            ha='center', va='center', color='black', fontsize=12)
    

# Draw vertical grid lines
for x in range(1, 3):
    axes[0].vlines(x, 0, 3, colors='black', linewidth=1)
    axes[1].vlines(x, 0, 3, colors='black', linewidth=1)

# Draw horizontal grid lines
for y in range(1, 3):
    axes[0].hlines(y, 0, 3, colors='black', linewidth=1)
    axes[1].hlines(y, 0, 3, colors='black', linewidth=1)


plt.savefig("Traffic lights for AE")