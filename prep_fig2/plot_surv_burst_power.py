import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import xarray as xr
import pandas as pd

print("reading survey...")
survey = pd.read_csv('/data/emfisis_burst/wip/rablack75/rablack75/CountSurvey/AllSurveyPowerA_Version20240515_withchorus_pos.csv')
# do the MLT fix
survey["MLT_rads"] = 2*np.pi*survey["MLT"]/24
# find the chorus in survey 
survey_chorus_hist = survey[(survey["chorus_pos"]==True) & (survey["Plasmapause"]=="Out") & (survey["Lstar"]>2.)]

print("reading burst...")
burst_chor = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/all_stats_mean_redidV12_fpefceADDED.nc')

print("plotting...")
plt.style.use('ggplot')
plt.figure(figsize=(6,5))
plt.title(r'$Integrated\ chorus\ power\ (f_{lhr} - f_{ce})$',fontsize=14)
sns.histplot(np.log10(survey_chorus_hist["Survey_integral"]),stat='probability',binwidth=0.1,color='darkgreen',label='Survey chorus events')
burst_pows = np.array(burst_chor["total_power"]).flatten()
burst_pows = burst_pows[burst_pows<10**5]
sns.histplot(np.log10(burst_pows),stat='probability',binwidth=0.1,color='lightgreen', label="Burst chorus events")
plt.xlabel(r'$Power\ (nT^{2})$',fontsize=14)
plt.ylabel(r'$Probability$',fontsize=14)
plt.xlim(-9,0)
plt.legend()

plt.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig2/attemptV1.png')