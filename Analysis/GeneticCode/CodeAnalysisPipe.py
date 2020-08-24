from Analysis.GeneticCode.CodonTransitions import create_codon_trans_matrix
from Data.config import SNP, General, CodeAnalysis
from lib.Utils import mkdirifnotexists
from os.path import join, basename
from Data.SNP.codons import get_codon_table
from pandas.io.pickle import read_pickle
import numpy as np
from Analysis.GeneticCode.CodonRisks import codon_risk
from glob import glob
from lib import Utils
from pandas.core.frame import DataFrame
from collections import defaultdict
import matplotlib.pyplot as plt
import Levenshtein
from scipy.stats.stats import mannwhitneyu
from scipy.stats.morestats import wilcoxon
from pandas.core.reshape.concat import concat

def million_codes():
    # This creates one million permutations of the genetic code
    aas, _, _ = get_codon_table()
    df = read_pickle(join(General.Basepath, 'All_4_60_mutation_codon_counts.df'))
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    for i in range(100):
        codon_risk(df, aas, 'All_{:02d}'.format(i), True, subdir='Million')
    compiled_f = join(CodeAnalysis.CodonsDir, 'Codon_risk_compiled.dat')
    ret = defaultdict(list)
    for i, fn in enumerate(glob(join(CodeAnalysis.CodonsDir, 'Million', '*.dat'))):
        ret_l = Utils.Load(fn)
        for var in ['n+_risk', 'c+_risk', 'o+_risk', 'hyd_risk', 'PR_risk']:
            ret[var].extend((ret_l[var] if i==0 else ret_l[var][1:]))
        print(i)
    Utils.Write(compiled_f, ret)
    return compiled_f

def plot_code_histograms(compiled_f, outdir):
    ret = Utils.Load(compiled_f)
    npr = np.array(ret['n+_risk'])
    cpr = np.array(ret['c+_risk'])
    hydr = np.array(ret['hyd_risk'])
    prr = np.array(ret['PR_risk'])
    for nm, riskarr, color in zip(['N_plus','C_plus','hyd', 'PR'], [npr, cpr, hydr, prr],
                                   ['#0d4c7c', '#151515', '#018571', '#660099']): 
        for nm1, riskarr1 in zip(['N_plus','C_plus','hyd', 'PR'], [npr, cpr, hydr, prr]):
            if nm==nm1:
                locarr = riskarr[1:]
                stan = riskarr[0]
            else:
                locarr = riskarr[1:][tuple([riskarr1[1:]<=riskarr1[0]])]
                stan = riskarr[0]
            _, ax = plt.subplots(1, figsize =(3.5,2.333), dpi=144)
            ax.hist(locarr, color=color, bins=100, density=True)
            ax.axvline(stan, color='yellow', lw=1)
            ax.axvline(stan, color='k', lw=0.6)
            print('{} given {} {} p={}'.format(nm, nm1, stan, sum(locarr<=stan)/len(locarr)))
            plt.savefig(join(outdir, 'Code_cost_million_hist_{}_{}.png'.format(nm, nm1)), dpi=144)
            plt.close('all')

def multi_organism_analyze():
    # This replicates the analysis presented in Fig. 4
    # Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5581930/
    codons_all = read_pickle('./resource/ModelOrganisms.df').set_index('Taxid')
    # Take only the organisms with more than 50K codons in the calculation
    codons_all = codons_all.loc[codons_all.iloc[:,11:].sum(1)>=50000]
    aas, _, _ = get_codon_table()
    # Create alternative codes for each organism and transition-transverskion rate
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    for taxid, row in codons_all.iterrows():
        codons = row[11:].astype(float)
        for titv in [0.2,0.25,0.333,0.5,0.667,1,1.5,2,3,4,5]:
            ti = (2*titv)/(1+2*titv)
            codon_risk(None, aas, 'Tax_{}_Rate_{}'.format(taxid, titv), all_mutations=False, 
                      external_counts=codons, external_titv=(ti,1-ti), subdir='MultiOrg')
    # Collate the results in one table
    proc_stats = {}
    for fnm in glob(join(CodeAnalysis.CodonsDir, 'MultiOrg/*.dat')):
        tax = float(basename(fnm).split('Tax_')[-1].split('_')[0])
        rate = float(basename(fnm).split('Rate_')[-1].split('_')[0])
        ret = Utils.Load(fnm)
        npr = np.array(ret['n+_risk'])
        cpr = np.array(ret['c+_risk'])
        proc_stats[(tax, rate)] = {'cpr_p':sum(cpr[1:]<=cpr[0])/10000., 
                                   'npr_p':sum(npr[1:]<=npr[0])/10000.,
                                   'ncpr_p': sum((cpr[1:]<=cpr[0]) & (npr[1:]<=npr[0])), 
                                   'cpr':cpr[0], 'npr':npr[0]}
    DataFrame(proc_stats).to_pickle(join(CodeAnalysis.CodonsDir, 'MultiOrg_rates.df'))

def issquare(cod):
    sq = True
    codons = [i for i in cod[(cod=='R') | (cod=='K') | (cod == 'Q')].index if i.endswith('A')]
    for i in codons:
        dsts = [Levenshtein.distance(i,j) for j in codons]
        sq = sq and sorted(dsts) == [0,1,1,2]
    return sq

def isdiag(cod):
    sq = True
    codons = [i for i in cod[(cod=='R') | (cod=='K') | (cod == 'Q')].index if i.endswith('A')]
    for i in codons:
        dsts = [Levenshtein.distance(i,j) for j in codons]
        sq = sq and sorted(dsts) == [0,2,2,2]
    return sq

def square_vs_diag(codon_permutation_f, outdir):
    # Replicates the analysis presented in Fig. 5B and creates plot
    ret = Utils.Load(codon_permutation_f)
    npr = np.array(ret['n+_risk'])
    codes = ret['code']
    squares = np.array([issquare(i) for i in codes[1:]])
    diags = np.array([isdiag(i) for i in codes[1:]])
    print ('n diags: {}'.format(sum(diags)))
    print ('n squares: {}'.format(sum(squares)))
    _, ax = plt.subplots(1, figsize =(4.7,5.2), dpi=144)
    grps = [npr[1:][squares], npr[1:][~squares&~diags], npr[1:][diags]]
    ax.boxplot(grps, showfliers=False,
                whis=(5,95), 
                flierprops={'color':'k', 'marker':'x', 'markersize':2}, 
                boxprops={'color':'k', 'lw':0.6}, capprops={'color':'k', 'lw':0.6},
                whiskerprops={'color':'k', 'lw':0.6},
                medianprops={'color':'', 'lw':1.2})
    ax.set_ylim(0.15,0.31)
    ax.set_yticks([0.15,0.2,0.25,0.3])
    print('squares vs all: {}'.format(mannwhitneyu(grps[0], grps[1])))
    print('squares vs diags: {}'.format(mannwhitneyu(grps[0], grps[2])))
    print('diags vs all: {}'.format(mannwhitneyu(grps[2], grps[1])))
    plt.savefig(join(outdir, 'Squares_diags.png'), dpi=144)

def codon_bias(outdir):
    # Data source: https://elifesciences.org/articles/41043
    df = read_pickle('./resource/Proc_strains_codons.df')
    thr = df.loc[['ACT','ACA','ACC','ACG']]
    thr = thr.truediv(thr.sum()).T
    ile = df.loc[['ATT','ATA','ATC']]
    ile = ile.truediv(ile.sum()).T
    print(wilcoxon(thr['ACT'], thr['ACA']))
    print(wilcoxon(thr['ACC'], thr['ACG']))
    print(wilcoxon(ile['ATT'], ile['ATA']))
    _, ax = plt.subplots(1, figsize =(8.2,5.2), dpi=144)
    bp = ax.violinplot(concat([thr, ile], axis = 1).T, positions=[1,2,3,4,6,7,8])
    for partname in ('cbars','cmins','cmaxes'):
        vp = bp[partname]
        vp.set_edgecolor('k')
        vp.set_linewidth(1)

    [m.set_color('#0d4c7c') for m in bp['bodies'][:4]]
    [m.set_color('#891919') for m in bp['bodies'][-3:]]
    ax.set_ylim(0,1.)
    ax.set_xticks(range(1,9))
    plt.savefig(join(outdir, 'Codon_usage.png'), dpi=144)

if __name__ == '__main__':
    create_codon_trans_matrix(SNP.OM_RGC.OutDir, mkdirifnotexists(join(SNP.OM_RGC.OutDir,'play')), 'All')
    
    # Replicates analysis presented in Fig. 2, 3
    # Replace outdir with a desired output directory
    compiled_f = million_codes()
    plot_code_histograms(compiled_f, outdir='.') 
    
    # Replicates analysis presented in Fig. 4 and saves it as a table
    multi_organism_analyze()
    
    # Replicates analysis in Fig. 5A and creates plot
    # Replace outdir with a desired output directory
    codon_bias(outdir='.')
    
    # Replicates analysis in Fig. 5B and creates plot
    # Replace codon_permutation_f with a file of codon permutations, e.g. one .dat file from
    # the 'Million' directory created above and outdir with a desired output directory
    square_vs_diag(codon_permutation_f=None, outdir=None) 