# This code creates the data required to make the bands of power plots (low frequency, lowerband and upperband)
# this time seperated by fpe/fce
# CURRENTLY this is using slightly incomplete data as the input whilst the full job is running


import xarray as xr
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.patches as patches


# other
import pandas as pd
import string
import roman

plt.style.use("ggplot")


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


def make_dial_plot(ax,histogram,colorbar,title,cmap):  
  
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
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)],fontsize=12)
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
    ax.set_xticklabels([str(i) if i % 6 == 0 else '' for i in range(0,24)],fontsize=12)
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

def regional_average(dist,region):
    ''' Takes the distributions of shape (240,70) (0.1MLT,0.1L),
        and averages between desired MLT sector '''
    
    if region == "dawn/post-mid":
        # between 1900 and 1200 (past midnight, so need to combine the two segments
        # for numpy array, as it doesn't wrap through 0)
        wrapped = np.concatenate((dist[190:], dist[:121]))
        average_power = np.nanmean(wrapped)
    
    return average_power


# import the chorus dataset
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/all_stats_mean_redidV12_fpefceADDED.nc')
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus270525.nc')
#chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus020625_fpefceadded.nc')
#chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

burstA = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/A/full_12_18_withTotalPowerAndSurv.nc')
burstB = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/ChorusinBurstandSurvey/B/full_12_19_withTotalPowerAndSurv.nc')
chorus_result = xr.concat([burstA, burstB], dim="x")

# name the different bands for loop
names = ["exohiss_band", "lowerband", "upperband"]

# Make a 4 x 3 plots (do: survey power, burst peak power, IQR, burst peak/survey)

# Create a figure

#chorus_result = chorus_result.where((chorus_result["isChorus"]),drop=True)
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24
plot_names = list(string.ascii_lowercase)
# Lowercase Roman numerals from 1 to 10
labels = [roman.toRoman(i).lower() for i in range(1, 19)]

fpe_fce_conds = [(chorus_result["fpe_fce"]<5),(chorus_result["fpe_fce"]>=5)&(chorus_result["fpe_fce"]<9),(chorus_result["fpe_fce"]>=9)]
MLAT_cond = (np.abs(chorus_result["MLAT"])<6.)
flags = ['lf','lb','ub']

survey_regional_average = np.zeros(3)
burst_regional_average = np.zeros(3)

for flag,name in zip(flags,names):
    # create condition figure
    fig = plt.figure(layout='constrained', figsize=(21, 10))
    subfigs = fig.subfigures(1, 6)
    for i,fpe_fce_cond in enumerate(fpe_fce_conds):

        axsLeft = subfigs[2*i].subplots(4, 1, subplot_kw={'projection': 'polar'}) 
        axsRight = subfigs[(2*i)+1].subplots(2, 1,) 

        
        print(f'fpe/fce condition {i} started...')
        # set conditions for desired AE range
        conditions = (chorus_result[f'{name}_peak']<1e4)&(chorus_result[f'{name}_peak'] !=0)&(chorus_result["Lstar"]>0.)&(chorus_result[f'{name}_surv_power'] !=0) & (chorus_result[f'{flag}_flag']==1) & MLAT_cond & fpe_fce_cond

        survey_vals = np.where(conditions,chorus_result[f'{name}_surv_power'][:].values, np.nan)
        iqr_vals =  np.where(conditions,chorus_result[f'{name}_iqr'][:].values, np.nan)
        #mad_vals = np.where(conditions,chorus_result[f'{name}_mad'][:].values, np.nan)
        peak_vals = np.where(conditions,chorus_result[f'{name}_peak'][:].values, np.nan)
        median_vals = np.where(conditions,chorus_result[f'{name}_med'][:].values, np.nan)
        ratio_vals = peak_vals/survey_vals
        prob_diff = (peak_vals - survey_vals)/(peak_vals+survey_vals)
        iqr_over_median = iqr_vals/median_vals
    
        #ratio_vals = np.where(conditions,ratio, np.nan)
    
        
        # make histograms
        peak_hist,mlt_l_dists = make_stat_histogram(chorus_result,peak_vals)
        peakratio_hist, mlt_l_dists = make_ratio_histogram(chorus_result,ratio_vals)
        surv_hist, mlt_l_dists = make_stat_histogram(chorus_result,survey_vals)
        iqr_hist,mlt_l_dists = make_stat_histogram(chorus_result,iqr_vals)
        #mad_hist,mlt_l_dists = make_stat_histogram(chorus_result,mad_vals)
        iqr_median_vals,mlt_l_dists = make_stat_histogram(chorus_result,iqr_over_median)

        print(f'{name} histograms made...')
        # Kaine's smoothing technique before plotting

        survey_smooth = np.zeros((24*10,70))
        ratio_smooth = np.zeros((24*10,70))
        peak_smooth = np.zeros((24*10,70))
        iqr_smooth = np.zeros((24*10,70))
        #mad_smooth= np.zeros((24*10,70))
        iqr_median_smooth = np.zeros((24*10,70))


        for k in range(24):
            survey_smooth[k*10:(k+1)*10,:] = surv_hist[k,:]

            peak_smooth[k*10:(k+1)*10,:] = peak_hist[k,:]

            ratio_smooth[k*10:(k+1)*10,:] = peakratio_hist[k,:]

            iqr_smooth[k*10:(k+1)*10,:] = iqr_hist[k,:]

            #mad_smooth[k*10:(k+1)*10,:] = mad_hist[k,:]

            iqr_median_smooth[k*10:(k+1)*10,:] = iqr_median_vals[k,:]
        
        survey_regional_average[i] = regional_average(survey_smooth,"dawn/post-mid")
        burst_regional_average[i] = regional_average(peak_smooth,"dawn/post-mid")


        # save all data to list
        datas = [survey_smooth,peak_smooth,ratio_smooth,iqr_median_smooth]

        # define normalisations
        peak_norm = mcolors.LogNorm(vmin=10**(-4), vmax=10**(-2))
        survey_norm = mcolors.LogNorm(vmin=10**(-4), vmax=10**(-2))
        ratio_norm = mcolors.LogNorm(vmin=1e1,vmax=1e4)
        iqr_norm = mcolors.LogNorm(vmin=10**(-6), vmax=10**(-2))
        iqr_median_norm = mcolors.LogNorm(vmin=10**(0), vmax=10**(2))
        # save all to list
        norms = [survey_norm,peak_norm,ratio_norm,iqr_median_norm]

        cmaps = ['plasma','plasma','seismic','plasma']

        # now do plot
        print(f'condition {i} plots begun...')
        
        for j,(data,norm,cmap) in enumerate(zip(datas,norms,cmaps)):
            
        
            title = labels[j] +'.'
            axis = axsLeft[j]
            axis.set_facecolor('whitesmoke')
            axis.set_title(title,loc='left',fontsize=18)

            if j == 2:
                dial_plot = make_dial_frac_plot(axis,data,norm,title)
                axis.set_xticklabels([str(i) if (i==18 or i==6) else '' for i in range(0,24)],fontsize=18)
            else:
                dial_plot = make_dial_plot(axis,data,norm,title,cmap)

            if (j == 0):
                
                axis.set_xticklabels([str(i) if (i==18 or i == 12 or i==6) else '' for i in range(0,24)],fontsize=18)
            
            
            elif (j == 3):
                axis.set_xticklabels([str(i) if (i==18 or i == 0 or i==6) else '' for i in range(0,24)],fontsize=18)

            
            else:
                axis.set_xticklabels([str(i) if (i==18 or i==6) else '' for i in range(0,24)],fontsize=18)

            if i == 0:
                box = axis.get_position()
                axis.set_position([box.x0 - 0.15, box.y0, box.width, box.height])
            else:
                box = axis.get_position()
                axis.set_position([box.x0 - 0.05, box.y0, box.width, box.height])

        
        survey_vals =[i for i in survey_vals if i != 0]
        survey_vals = [i for i in survey_vals if i <1e4]
        log_bins = np.logspace(np.log10(np.nanmin(survey_vals)), np.log10(np.nanmax(survey_vals)), 50)
        print(log_bins)
        hist_survey = make_log_histogram(survey_vals,log_bins)
        hist_survey_bins = log_bins

        sns.histplot(data=survey_vals,bins=log_bins,color='lightblue',label="Survey power", kde=False,stat='probability',fill=False,element="step",ax=axsRight[0])

        peak_vals =[i for i in peak_vals if i != 0]
        peak_vals = [i for i in peak_vals if i <1e4]
        log_bins = np.logspace(np.log10(np.nanmin(peak_vals)), np.log10(np.nanmax(peak_vals)), 50)
        hist_peak = make_log_histogram(peak_vals,log_bins)
        hist_peak_bins = log_bins

        sns.histplot(data=peak_vals,bins=log_bins,color='darkblue',label="Peak power", kde=False,stat='probability',fill=False,element="step",ax=axsRight[0])

        title = labels[j+1] +'.'
        axsRight[0].axvline(0.0225,color='red',linestyle='dashed',label='Nonlinear threshold')
        axsRight[0].set_title(title,loc='left',fontsize=18)
        axsRight[0].set_xscale('log')
        sns.set(font_scale=0.8)
        axsRight[0].set_ylabel(r'Probabillity',fontsize=18)
        axsRight[0].set_xlabel(r'Power $(nT^2)$',fontsize=18)
        if (2*i)+1 == 1:
            axsRight[0].legend(loc="upper left")
        axsRight[0].set_xlim(1e-10,1e2)
        axsRight[0].set_facecolor('whitesmoke')


        title = labels[j+2] +'.'

        # now do the iqr probabillity
        # Now do the right hand side 
        #iqr_vals = [i for i in iqr_over_median if i != 0]
        #iqr_vals = [i for i in iqr_over_median if i <1e10]
        
        #log_bins = np.logspace(np.log10(np.nanmin(iqr_vals)), np.log10(np.nanmax(iqr_vals)), 50)
        #print(log_bins)
        #hist_iqr = make_log_histogram(iqr_vals,log_bins)
        #hist_iqr_bins = log_bins
        #sns.histplot(data=iqr_vals,bins=log_bins,color='orange',label="Burst power IQR/median", kde=False,stat='probability',element="step",fill=False,ax=axsRight[1])

        ratio_vals = [i for i in ratio_vals if i >10**(0)]
        ratio_vals = [i for i in ratio_vals if i <1e10]
        
        log_bins = np.logspace(np.log10(np.nanmin(ratio_vals)), np.log10(np.nanmax(ratio_vals)), 50)
        print(log_bins)
        sns.histplot(data=ratio_vals,bins=log_bins,color='green',label="Burst peak power/survey power", kde=False,stat='probability',element="step",fill=False,ax=axsRight[1])

        axsRight[1].set_title(title,loc='left',fontsize=18)
        axsRight[1].set_xscale('log')
        sns.set(font_scale=0.8)
        axsRight[1].set_ylabel(r'Probabillity',fontsize=18)
        axsRight[1].set_xlabel(r'Power ratio',fontsize=18)
        if (2*i)+1 == 1:
            axsRight[1].legend(loc="upper left")
        axsRight[1].set_facecolor('whitesmoke')

        axsRight[1].set_xlim(10**(0),10**(4))

        if i!=0:
            box = axsRight[0].get_position()
            axsRight[0].set_position([box.x0 - 0.05, box.y0, box.width, box.height])
            box = axsRight[1].get_position()
            axsRight[1].set_position([box.x0 - 0.05, box.y0, box.width, box.height])
        else:
            box = axsRight[0].get_position()
            axsRight[0].set_position([box.x0+0.05, box.y0, box.width, box.height])
            box = axsRight[1].get_position()
            axsRight[1].set_position([box.x0+0.05, box.y0, box.width, box.height])


    subfigs[0].suptitle(r'$a)\ f_{pe}/f_{ce}<5$', x=0.0, ha='left',fontsize=18)
    subfigs[2].suptitle(r'$b)\ 5\leq f_{pe}/f_{ce}<9$',x=0.0, ha='left', fontsize=18)
    subfigs[4].suptitle(r'$c)\ f_{pe}/f_{ce}\geq 9$', x=0.0, ha='left',fontsize=18)

    # add the colorbar to the leftside
    sfig = subfigs[0]
    # first, for power:
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=peak_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.7, 0.6, 0.05, 0.2])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1)
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    cbar.set_label(r'$nT^2$', fontsize=18) 

    # Add a separate axis for the color bar of the ratio
    sm = plt.cm.ScalarMappable(cmap='Greens', norm=ratio_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax2 = sfig.add_axes([0.7, 0.31, 0.05, 0.15])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax2.spines.values():
        spine.set_visible(False)
    cbar_ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar2=sfig.colorbar(sm, ax=cbar_ax2, orientation='vertical', fraction=0.5, pad=0.05,ticks=[1e1, 1e2, 1e3])
    cbar2.ax.tick_params(labelsize=14)  # Set colorbar tick label font size


    # Add colorbar for the iqr
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=iqr_median_norm)
    sm.set_array([])  # Required for the colorbar to work properly
    cbar_ax = sfig.add_axes([0.7, 0.11, 0.05, 0.15])  # [left, bottom, width, height]
    # Make the axes invisible
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    cbar_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add the color bar to the figure
    cbar=sfig.colorbar(sm, ax=cbar_ax, orientation='vertical', fraction=0.5, pad=0.1,ticks=[1e0, 1e2])
    cbar.set_ticklabels([r'$10^0$',r'$10^2$'])
    cbar.ax.tick_params(labelsize=14)  # Set colorbar tick label font size
    #cbar.set_label(r'$nT^2$', fontsize=14)

    print(f'{name} done')
            
    fig.savefig(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig5/all_latitudes/{name}_allLats_withratio.png')

    averages = [survey_regional_average,burst_regional_average]

    np.savetxt(f'/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/data_fig4/fpefce/combinedAB/LowMLATAverages_{name}.txt', averages, fmt="%.7f") 