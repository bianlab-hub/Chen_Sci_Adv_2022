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


mm10= bioframe.fetch_chromsizes('mm10')
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


#Make saddle plot
from cooltools import saddle
condx, condy = 'MC', 'FL'

QUANTILE_BINNING = True

binedges = {}
digitized = {}
hist = {}
sums = {}
counts = {}
saddledata = {}

#make a plot of 2 columns
gs = plt.GridSpec(nrows=1, ncols=2)
fig = plt.figure(figsize=(12, 6))
histbins = 30


for i, cond in enumerate([condx, condy]):
    exp = pd.read_table(f'./{long_names[cond]}.{binsize//1000}kb.expected.cis.tsv')
    eig = pd.read_table(f'./{long_names[cond]}.{binsize//1000}kb.eigs.cis.vecs.txt')
    # Determine how to bin the range of the E1 signal
    if QUANTILE_BINNING:
        q_binedges = np.linspace(0, 1, histbins)
        binedges[cond] = saddle.quantile(eig['E1'], q_binedges)
    else:
        qlo, qhi = saddle.quantile(eig['E1'], [0.02, 0.98])  # trim outliers
        binedges[cond] = np.linspace(qlo, qhi, histbins)

    digitized[cond], hist[cond] = saddle.digitize_track(binedges[cond], track=(eig, 'E1'))
    getmatrix = saddle.make_cis_obsexp_fetcher(clrs[cond], (exp, 'balanced.avg'))
    sums[cond], counts[cond] = saddle.make_saddle(getmatrix, binedges[cond],(digitized[cond], 'E1.d'),contact_type='cis')
    saddledata[cond] = sums[cond] / counts[cond]
    # Make the saddle plot
    g = saddle.saddleplot(
        q_binedges if QUANTILE_BINNING else binedges[cond],
        hist[cond], 
        np.log10(saddledata[cond]), 
        color=colors[cond],
        heatmap_kws={'vmin': -0.5, 'vmax': 0.5}, 
        fig=fig, subplot_spec=gs[i])

plt.savefig(condx+ '_' +condy + "_saddle_plot.pdf", dpi=300, format='pdf')

strength = {
    cond: saddle.saddle_strength(sums[cond], counts[cond]) 
        for cond in [condx, condy]
}
tmp = pd.DataFrame.from_dict(strength) 
tmp.to_csv(condx+ '_' +condy + "_saddle_strength.tsv", sep='\t', index=False)

gs = plt.GridSpec(nrows=1, ncols=2)
plt.figure(figsize=(14, 6))

plt.subplot(gs[0])
x = np.arange(histbins + 2)
for cond in [condx, condy]:
    plt.step(x[:-1], strength[cond], where='pre', color=colors[cond], label=long_names[cond])

plt.legend()
plt.xlabel('extent')
plt.ylabel('(AA + BB) / (AB + BA)')
plt.title('saddle strength profile')
plt.axhline(0, c='grey', ls='--', lw=1)
plt.xlim(0, len(x)//2)

plt.subplot(gs[1])
plt.step(x[:-1], strength[condy] / strength[condx], where='pre', c='k')
plt.axhline(1, c='grey', ls='--', lw=1)
plt.xlim(0, len(x)//2)
plt.xlabel('extent')
plt.ylabel('enrichment')
plt.title(condy + '/' + condx)
plt.savefig(condx+ '_' +condy + "_saddle_strength_profile.pdf", dpi=300, format='pdf')
