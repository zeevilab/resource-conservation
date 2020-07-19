import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import mannwhitneyu
from pandas.core.series import Series
from statsmodels.stats import multitest as smm
from pandas.io.pickle import read_pickle
from Analysis.GeneData.codons import get_codon_table
from RawData.config import Biodata, General, SNP
from pandas.core.frame import DataFrame
from os.path import basename, join, exists
from graph.plotutils import PdfIterator
from pandas.core.reshape.concat import concat
from itertools import product
from pandas.io.parsers import read_csv

def test_sig(loN, hiN, picout):
    num_samps = min(loN.shape[0], hiN.shape[0])//2
    loN = loN.loc[np.random.permutation(loN.index)]
    hiN = hiN.loc[np.random.permutation(hiN.index)]
    ret = {'lo':{}, 'hi':{}, 'lohi':{}}
    for col in loN:
        try:
            ret['lo'][col] = mannwhitneyu(loN[:num_samps][col], loN[num_samps:][col])
            ret['hi'][col] = mannwhitneyu(hiN[:num_samps][col], hiN[num_samps:][col])
        except ValueError:
            pass
    loN = loN.loc[np.random.permutation(loN.index)]
    hiN = hiN.loc[np.random.permutation(hiN.index)]
    for col in loN:
        try:
            ret['lohi'][col] = mannwhitneyu(loN[:num_samps][col], hiN[num_samps:][col])
        except ValueError:
            pass
    _, axes = plt.subplots(1,3, figsize=(12,5))  
    for i, rk in enumerate(ret.keys()):
        ldf = Series(ret[rk]).apply(Series).rename(columns={0:'U',1:'p'})
        ldf['q'] = smm.multipletests(ldf.p,alpha=0.05,method='fdr_bh')[1]
        minx = int(np.log10(ldf.q.min()))-1
        ldf['q'].hist(bins = np.logspace(minx,0,50), ax=axes[i])
        axes[i].set_xscale('log', basex=10)
        axes[i].set_title(rk)
    plt.savefig(picout)
        
        
def descriptive_codons(fname, norm_to_syn=False, mutation_flux=True):
    df = read_pickle(fname.replace('s.df','_counts.df')) if (norm_to_syn or \
                                                             (mutation_flux and not 'CLR' in fname))\
             else read_pickle(fname)
    cnts = read_pickle(fname.replace('mutation_codons','nucleotide_count').replace('CLR',''))
    aas, _, _ = get_codon_table()
    md = read_pickle(Biodata.United.SampleMeasurementsDF)
    df = (df.loc[:,cnts>1e6] if 'All' not in fname else df.loc[:,cnts>1e7])
    if norm_to_syn:
        def normfunc(ldf):
            ldf = ldf.drop((ldf.name, ldf.name))
            syn = [i for i in ldf.index if aas[i[0]] == aas[i[1]]]
            if len(syn) == 0:
                return
            ldf = ldf.truediv(ldf.loc[syn].mean()).drop(syn)
            ldf.index = ldf.index.get_level_values(1)
            return ldf
        df = df.groupby(level=0).apply(normfunc)
    if mutation_flux:
        if not norm_to_syn and not 'CLR' in fname:
            #CLR
            df = df.add(1).apply(np.log10)
            df = df.subtract(df.mean())
        ldf = {}
        for nm,row in df.iterrows():
            if (nm[1],nm[0]) in df.index:
                if norm_to_syn:
                    ldf[nm] = np.log(row.truediv(df.loc[(nm[1],nm[0])]))
                else:
                    ldf[nm] = row.subtract(df.loc[(nm[1],nm[0])])
        ldf = DataFrame(ldf)
#         df = ldf.truediv(cnts, axis=0).T
        df = ldf.T
    dfx = df.T.join(md[['Nitrate_umolkg-1', 'Depth_m']]).dropna(subset=['Nitrate_umolkg-1', 'Depth_m'])\
            .sort_values(['Nitrate_umolkg-1', 'Depth_m']).drop(['Nitrate_umolkg-1', 'Depth_m'], axis=1)
    test_sig(dfx[:80], dfx[-80:], '/ru-auth/local/home/dzeevi/ndiff/codonmuts_' + basename(fname).split('_')[0] + \
                         ('_syn' if norm_to_syn else '')\
                      + ('_logodds' if mutation_flux else '') + '_sigs.png')
    loN = dfx[:40]
    hiN = dfx[-40:]
    ret = {}
    for c in loN.columns:
        try:
            mwu = mannwhitneyu(loN[c].dropna(), hiN[c].dropna())
        except:
            continue
        med_diff = hiN[c].median() - loN[c].median()
        ret[c] = {'U':mwu[0], 'Diff':med_diff, 'p':mwu[1]}
    ser = DataFrame(ret).T
    ser['q'] = smm.multipletests(ser.p,alpha=0.05,method='fdr_bh')[1]
    ser['aa_start'] = [aas[i] for i in ser.index.get_level_values(0)]
    ser['aa_end'] = [aas[i] for i in ser.index.get_level_values(1)]
    from pandas import read_csv
    aa_props = read_csv('/ru-auth/local/home/dzeevi/tmp/aa_s.csv').set_index('Code')
    del read_csv
    ser = ser.join(aa_props.rename(columns = {c:'start_'+c for c in aa_props}), on='aa_start')\
            .join(aa_props.rename(columns = {c:'end_'+c for c in aa_props}), on='aa_end')
    for col in ser.columns:
        if col.startswith('start_'):
            ser[col.replace('start','diff')] = ser[col.replace('start','end')]-ser[col]
    ser.to_csv('/ru-auth/local/home/dzeevi/ndiff/codonmuts_' + basename(fname).split('_')[0] + \
                         ('_syn' if norm_to_syn else '')\
                      + ('_logodds' if mutation_flux else '') + '.csv')
    pdf = PdfIterator('/ru-auth/local/home/dzeevi/ndiff/codonmuts_' + basename(fname).split('_')[0] + \
                         ('_syn' if norm_to_syn else '')\
                      + ('_logodds' if mutation_flux else ''),  
                      size=(9,16), colfigs=1, rowfigs=3)
    longcnct = []
    for _, grp in ser.groupby(level=0):
        tocnct=[]
        for nm, row in grp.iterrows():
            cnct = tocnct if nm[0] != nm[1] else longcnct
            cnct.append(loN[nm].to_frame('loN_{}({})_to_{}({}) q={:.1e}'.format(nm[0],row['aa_start'],
                                                                                nm[1], row['aa_end'], row['q'])))
            cnct.append(hiN[nm].to_frame('hiN_{}({})_to_{}({}) q={:.1e}'.format(nm[0],row['aa_start'],
                                                                                nm[1], row['aa_end'], row['q'])))
        ax = next(pdf.iter)
        concat(tocnct, axis=1, sort=False).boxplot(rot=90,ax=ax)
        if not mutation_flux:
            ax.set_yscale('log', basey=10)
    if not norm_to_syn and len(longcnct)>0:
        ldf = concat(longcnct, axis=1, sort=False)
        for i in range(0,ldf.shape[1], 20):
            ax = next(pdf.iter)
            ldf.iloc[:,i:i+20].boxplot(rot=90,ax=ax)
            if 'CLR' not in fname and not mutation_flux:
                ax.set_yscale('log', basey=10)
    pdf.Terminate()
#         print('jackie')
#     ser.q.hist(bins= np.logspace(-17,0,35))
#     plt.xscale('log', basex=10)
#     plt.show()
    print('jackie')
    
def calc_codon_costs(out_f=None, force_rerun=False):
    if exists(out_f) and not force_rerun:
        return read_pickle(out_f)
    aas, _, _ = get_codon_table()
    ret = {}
    for cod_1, aa1 in aas.to_dict().items():
        for cod_2, aa2 in aas.to_dict().items():
            ret[(cod_1, cod_2)] = {'aa_s':aa1, 'aa_e':aa2}
    aa_props = read_csv('/ru-auth/local/home/dzeevi/tmp/aa_NCHP.csv', index_col=0)
    ret = DataFrame(ret).T
    ret.index.names = ['Codon_s','Codon_e']
    codon_props = ret.join(aa_props, on='aa_s').join(aa_props, on='aa_e', lsuffix='_s', rsuffix='_e')
    for v in ['C','N','hyd','PR']:
        codon_props[v+'_d'] = codon_props[v+'_e']-codon_props[v+'_s']
        codon_props[v+'_abs_d'] = abs(codon_props[v+'_d'])
    if out_f is not None:
        codon_props.to_pickle(out_f)
    return codon_props
    
if __name__=='__main__':
    calc_codon_costs(SNP.CodonProps, True)
#     descriptive_codons(join(General.Basepath, 'AllCLR_4_60_mutation_codons.df'), False, True)
#     for a,b in product([False, True],[False,True]):
#         descriptive_codons(join(General.Basepath, 'Top100_4_60_mutation_codons.df'), a, b)