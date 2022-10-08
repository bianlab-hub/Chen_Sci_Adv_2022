import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocess as mp
import bioframe
import cooler
import itertools
import click
import cooltools
import cooltools.eigdecomp
import cooltools.expected
import cooltools.saddle
from dask.distributed import Client, LocalCluster
from scipy.linalg import toeplitz
import scipy.ndimage



# regions in the format: reg1=('chr1', 40000000, 45000000)

reg1=('chr11', 111000000, 113500000)
chr=reg1[0]
sample='FL_combined_test'
#set fixed min max for observed matrix
obsmin=-4
obsmax=-4.8

oemin=-1
oemax=1


coolfile='/data05/CQM/CQM_20200401_HiC_combined/Forelimb_cool/E145_Forelimb_combined_10kb.cool'
cisexpfile='/data05/CQM/CQM_20200401_HiC_combined/Forelimb_cool/E145_Forelimb_combined_cis_expected.tsv'
c = cooler.Cooler(coolfile)
cis_expected=pd.read_csv(cisexpfile, sep='\t')
cis_expected_2 = {k: x.values for k, x in cis_expected.groupby('chrom')['balanced.avg']}
obs_mat = c.matrix().fetch(reg1)
exp_mat = toeplitz(cis_expected_2[chr][:obs_mat.shape[0]])

#normalize with total interactions in the region
#Change1:nansum(obs_mat)>>np.nansum(obs_mat)
obs_mat_norm = obs_mat/(np.nansum(obs_mat))

#plot obs heatmap 
#Change2:obs_mat >> obs_mat_norm   
mat=np.log10(obs_mat_norm)
row_chrom=chr
col_chrom=chr
scale='log10'
out=sample + '_' + reg1[0]+'_'+ str(reg1[1]) + '_' + str(reg1[2]) + '_obs_10kb'+ '.pdf'
dpi= 300
colormap='YlOrRd'
#obsmin=np.log10(np.nanquantile(obs_mat_norm,0.05))
#obsmax=np.log10(np.nanquantile(obs_mat_norm,0.95))

plt.figure(figsize=(10,10))
plt.gcf().canvas.set_window_title("Contact matrix".format())
plt.title(sample + " obs heatmap")
plt.imshow(mat, interpolation="none",vmin=obsmin,vmax=obsmax, cmap=colormap)
plt.ylabel("{} coordinate".format(row_chrom))
plt.xlabel("{} coordinate".format(col_chrom))
#label axis with Mbs
BINSIZE=10000
ticks_pixels = np.linspace(0, (reg1[2]-reg1[1])//BINSIZE,((reg1[2]-reg1[1])//1000000+1))
ticks_mbp = ((ticks_pixels*BINSIZE + reg1[1])//1000000).astype(int)
plt.xticks(ticks_pixels, ticks_mbp)
plt.yticks(ticks_pixels, ticks_mbp) 
## Label PRS/Proximal TBC/Sox9 TSS region with plt.axvline 
# PRS: chr11	111555526	111833044	
# Proximal TBC chr11 112437850	112722603
# Sox9 TSS chr11 112782209
plt.axvline(x=55.5526,ls="dashed",c="red")
plt.axvline(x=83.3044,ls="dashed",c="red")
plt.axvline(x=143.785,ls="dashed",c="blue")
plt.axvline(x=172.2603,ls="dashed",c="blue")
plt.axvline(x=178.22093,ls="dashed",c="grey")
cb = plt.colorbar()
cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
    }[scale])
plt.savefig(out, dpi=dpi, format='pdf')
plt.close()

#plot oe heatmap
#Change3:obs_mat >> obs_mat_norm 
mat=np.log2(obs_mat_norm/exp_mat)
row_chrom=chr
col_chrom=chr
scale='log2'
out=sample + '_' + reg1[0]+'_'+ str(reg1[1]) + '_' + str(reg1[2]) + '_oe_10kb'+ '.pdf'
dpi= 300
colormap='bwr'

1
plt.figure(figsize=(10,10))
plt.gcf().canvas.set_window_title("Contact matrix".format())
plt.title(sample + " obs/exp heatmap")
plt.imshow(mat, interpolation="none",vmin=oemin,vmax=oemax, cmap=colormap)
plt.ylabel("{} coordinate".format(row_chrom))
plt.xlabel("{} coordinate".format(col_chrom))
ticks_pixels = np.linspace(0, (reg1[2]-reg1[1])//BINSIZE,((reg1[2]-reg1[1])//1000000+1))
ticks_mbp = ((ticks_pixels*BINSIZE + reg1[1])//1000000).astype(int)
plt.xticks(ticks_pixels, ticks_mbp)
plt.yticks(ticks_pixels, ticks_mbp)
cb = plt.colorbar()
cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
    }[scale])
plt.savefig(out, dpi=dpi, format='pdf')
plt.close()

