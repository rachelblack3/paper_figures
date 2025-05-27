import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import seaborn as sns
import pandas as pd

# functions for making the power dial plots
plt.style.use("ggplot")

chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus_burst260525.nc')
# save total power to chorus_result
burst_pow = chorus_result["mean_burst_int"].where((np.abs(chorus_result["MLAT"])>=6.))
print("burst read...")

survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerA_Version20240515_withchorus_pos.csv')
survey["MLT_rads"] = 2*np.pi*survey["MLT"]/24
# find the chorus in survey 
survey_chorus = survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>3.) & (np.abs(survey["MLAT"])>=6.)]
survey_pow = survey_chorus["Survey_integral"]
# taking a 1 L* and 2 MLT bin:
chorus_result["MLT_rads"] = 2*np.pi*chorus_result["MLT"]/24

# do the plot
fig, axes = plt.subplots(6, 3, figsize=(4*3, 7*3),sharey=True)

for j,MLT in enumerate([0,2,4,6,8,10]):
    
    print(f'The MLT range for the following plots are {MLT} to {MLT+2}')

    for i,Lstar in enumerate([3,4,5]):
        ax = axes[j,i]
        
            
        print(f'The L* range for the following plots are {Lstar} to {Lstar+1}')
        surv_scen1 = survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>Lstar) & (survey["Lstar"]<=Lstar+1) & (survey["MLT"]>MLT) & (survey["MLT"]<=MLT+2)&(np.abs(survey["MLAT"])>=6.)]
        surv_pow_scen1 = surv_scen1["Survey_integral"]

        burst_pow_scen1 = chorus_result["mean_burst_int"].where((np.abs(chorus_result["MLAT"])<6.)& (chorus_result["Lstar"]>Lstar) & (chorus_result["Lstar"]<=Lstar+1) & (chorus_result["MLT"]>MLT) & (chorus_result["MLT"]<=MLT+2))

        log_bins = np.linspace(-8, 0, 40)

        sns.histplot(np.log10(surv_pow_scen1),stat='probability',bins=log_bins,color='darkgreen', label="Survey",ax=ax)
        sns.histplot(np.log10(burst_pow_scen1),stat='probability',bins=log_bins,color='lightgreen', label="Burst",ax=ax)

        if j ==5:
            ax.set_xlabel(r'$Chorus\ average\ power\ (nT^2)$',fontsize=14)
        else:
            ax.set(xlabel=None)
            ax.set_xticklabels(labels=[])
        #if i == 0:
        #    ax.set_ylabel(r'$Probabillity$',fontsize=14)
        #else:
        #    ax.set(ylabel=None)
        #    ax.set_yticks(ticks=[],labels=[])
            
        sns.set(font_scale=0.8)
        ax.set_xlim(-8,0)

axes[0,0].legend()

for j in range(6):
    MLT = 0
    
    axes[j,0].set_ylabel(r'$Probabillity$',fontsize=14)
    axes[j,0].text(-0.65, 0.5, f'{(MLT+j*2)}-{MLT+(j+1)*2} MLT', transform=axes[j,0].transAxes,
        ha='left', fontsize=12)


    axes[j,1].set(ylabel=None)
    axes[j,2].set(ylabel=None)
    axes[j,1].tick_params(labelleft=False)
    axes[j,2].tick_params(labelleft=False)  

for i in range(3):
    L_star = i+3
    
    axes[0,i].text(0.5, 1.05, f'{L_star}-{L_star+1} L*', transform=axes[0,i].transAxes,
        ha='center', fontsize=12)   
    
# Add extra space on the left
fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05,
                    wspace=0.1, hspace=0.4)

fig.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/HighLatSplitDist_Part1.png')