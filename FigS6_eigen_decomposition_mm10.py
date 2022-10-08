from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')

import multiprocess as mp
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler



mm10 = bioframe.fetch_chromsizes('mm10')
chromsizes = bioframe.fetch_chromsizes('mm10')
chromosomes = list(chromsizes.index)


conditions = ['MC', 'FL']
binsize = 10000

cooler_paths = {    
    'MC' : f'/data05/CQM/CQM_20200401_HiC_combined/E145_MC_combined_10kb.cool',
    'FL' : f'/data05/CQM/CQM_20200401_HiC_combined/E145_Forelimb_combined_10kb.cool',
}
long_names = {
    'MC': 'E14.5 MC',
    'FL' : 'E14.5 FL',
}
pal = sns.color_palette('colorblind')
colors = {
    'MC': pal[0],
    'FL' :  pal[1],
}

clrs = {
    cond: cooler.Cooler(cooler_paths[cond]) for cond in conditions
}



bins = cooler.binnify(mm10, binsize)
fasta_records = bioframe.load_fasta('/data05/genomes/mm10_20chr.fa')
bins['GC'] = bioframe.tools.frac_gc(bins, fasta_records)
bins.head()

bins2 = cooler.binnify(clrs['MC'].chromsizes, binsize)
bins2['GC'] = bioframe.tools.frac_gc(bins2, fasta_records)
bins2.head()

_=plt.hist(bins['GC'].dropna(), range=(0.2, 0.6), bins=100)
plt.xlabel('GC content')
plt.title(f'mm10, {binsize//1000}kb bins')
plt.savefig("mm10_GC_content.pdf", dpi=300, format='pdf')

from cooltools.eigdecomp import cooler_cis_eig

lam = {}
eigs = {}

for cond in conditions:
    lam[cond], eigs[cond] = cooler_cis_eig(
        clrs[cond], 
        bins,
        n_eigs=3, 
        phasing_track_col='GC', 
        sort_metric='var_explained')
    
    # Save text files
    lam[cond].to_csv(f'./{long_names[cond]}.{binsize//1000}kb.eigs.cis.lam.txt', sep='\t')
    eigs[cond].to_csv(f'./{long_names[cond]}.{binsize//1000}kb.eigs.cis.vecs.txt', sep='\t', index=False)
    
    # Save bigwig track
    bioframe.to_bigwig(eigs[cond], mm10, f'./{long_names[cond]}.{binsize//1000}kb.eigs.cis.vecs.E1.bw', 'E1')

#correlation between E1 values
from scipy.stats import rankdata

gs = plt.GridSpec(nrows=1, ncols=2)
plt.figure(figsize=(16, 6))
condx, condy = 'MC', 'FL'

plt.subplot(gs[0])
lo, hi = -2 , 2
plt.hexbin(
    eigs[condx]['E1'],
    eigs[condy]['E1'],
    vmax=50,
)
plt.xlabel('E1 ' + long_names[condx])
plt.ylabel('E1 ' + long_names[condy])
plt.gca().set_aspect(1)
plt.xlim(lo, hi)
plt.ylim(lo, hi)
plt.axvline(0, c='b', lw=0.5, ls='--')
plt.axhline(0, c='b', lw=0.5, ls='--')
plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
plt.colorbar(shrink=0.6)
plt.savefig("MC_FL_E1_correlation_plot.pdf", dpi=300, format='pdf')


plt.subplot(gs[1])
mask = eigs[condx]['E1'].notnull() & eigs[condy]['E1'].notnull() 
vx = eigs[condx]['E1'].loc[mask].values
vy = eigs[condy]['E1'].loc[mask].values
lo, hi = 0 , len(vx)

plt.hexbin(
    rankdata(vx),
    rankdata(vy),
    vmax=20,
)
plt.xlabel('E1 rank ' + long_names[condx])
plt.ylabel('E1 rank ' + long_names[condy])
plt.gca().set_aspect(1)
plt.xlim(lo, hi)
plt.ylim(lo, hi)
plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
plt.colorbar(shrink=0.6)
plt.savefig("MC_FL_E1_rank_correlation_plot.pdf", dpi=300, format='pdf')
