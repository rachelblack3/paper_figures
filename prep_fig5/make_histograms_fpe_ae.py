import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


plt.style.use("ggplot")

chorus_result = xr.open_dataset('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/burst_analysed/just_chorus020625_fpefceadded.nc')

fpe_fce_vals = chorus_result["fpe_fce"]
AE = chorus_result["AE"]

fig,axs= plt.subplots(1,2, figsize=(10,6),sharey=True,gridspec_kw={'wspace': 0.1})

axs[0].hist(AE[AE<=1e4],bins=25)
axs[0].set_xlabel(r'$AE$',fontsize=14)
axs[0].set_ylabel(r'$Count$',fontsize=14)
axs[0].set_xlim(0,2000)

axs[1].hist(fpe_fce_vals[fpe_fce_vals<=25.],bins=25)
axs[1].set_xlabel(r'$f_{pe}/f_{ce}$',fontsize=14)

fig.savefig('/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plots_fig5/fpece_ae_dists.png')