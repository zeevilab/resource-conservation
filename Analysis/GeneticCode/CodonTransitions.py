from glob import glob
from os.path import join, basename, exists
from pandas.io.pickle import read_pickle
from pandas.core.series import Series
from pandas.core.reshape.concat import concat
from pandas.core.indexes.multi import MultiIndex
from itertools import product
from lib.Utils import chdirmkifnotexist
from Data.config import General
from Data.SNP.codons import get_codon_table

def get_relevant_codonpairs():
    lst = []
    for i,j,k in product('ACGT',repeat=3):
        oldc = i+j+k
        for l in range(3):
            for m in 'ACGT':
                if m != oldc[l]:
                    newc =oldc[:l]+m+oldc[l+1:]
                    lst.append((oldc,newc))
    return lst

def do_one_unite(fname, grpnm, tmpdir):
    cod = read_pickle(fname)
    fcod = read_pickle(fname.replace('codons','percmut'))
    allcodons = get_relevant_codonpairs()
    def apply_counts(ldf):
        nm = ldf.index.get_level_values(0)[0]
        ldf.index = MultiIndex.from_tuples(ldf.index.get_level_values(1))
        ldf = ldf.dropna(how='all', axis=1)
        codldf = Series({i:(cod.loc[nm, i[0]] if i[0] in cod.loc[nm] else 0) for i in allcodons})
        ldf = codldf.to_frame('cod').join(ldf).drop('cod', axis=1)
        for col in ldf.columns:
            ldf[col + '_s'] = codldf
        return ldf.fillna(0)
    fcod = fcod.groupby(level=0).apply(apply_counts)
    fcod.groupby(level=[1,2]).sum().to_pickle(join(tmpdir, grpnm+basename(fname).replace('codons','sums')))
    
def create_codon_trans_matrix(dirsname, tmpdir, grpnm):
    chdirmkifnotexist(join(tmpdir, 'tmpmq'))
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    for fname in glob(join(dirsname, '*codons*')):
        if exists(join(tmpdir, grpnm+basename(fname).replace('codons','sums'))):
            continue
        do_one_unite(fname, grpnm, tmpdir)
    df = None
    for f in glob(join(tmpdir, grpnm+'*sums.df')):
        print(basename(f))
        if df is None:
            df = read_pickle(f)
        else:
            df = concat([df,read_pickle(f)], sort=False).groupby(level=[0,1]).sum()
    df.columns = MultiIndex.from_tuples([('per',i) if '_s' not in i \
                                         else ('sums',i.replace('_s','')) for i in df.columns])
    sms = df.sums.sum().truediv(3)
    def applyfunc(ldf):
        ret = ldf.per.truediv(ldf.sums).loc[ldf.name].iloc[:,0]
        ret[ldf.name] = 1-ret.sum()
        return ret
    def applyfunc_counts(ldf):
        adds = ldf.sums.iloc[0][0]
        ret = ldf.per.loc[ldf.name].iloc[:,0]
        ret[ldf.name] = adds - ret.sum()
        return ret
    df_counts = df.groupby(level=1, axis = 1).apply(lambda x: x.groupby(level=0).apply(applyfunc_counts))
    df_norm = df.groupby(level=1, axis = 1).apply(lambda x: x.groupby(level=0).apply(applyfunc))
    aas, _, _ = get_codon_table()
    df_aas = df_counts
    df_aas['aa_start'] = [aas[a] for a in df_counts.index.get_level_values(0)]
    df_aas['aa_end'] = [aas[a2] + ('.' if a1!=a2 and aas[a1]==aas[a2] else '') \
                        for a1, a2 in \
                        zip(df_counts.index.get_level_values(0),df_counts.index.get_level_values(1))]
    df_aas = df_aas.groupby(['aa_start','aa_end']).sum()
    df_norm.to_pickle(join(General.Basepath, grpnm + '_4_60_mutation_codons.df'))
    df_counts.to_pickle(join(General.Basepath, grpnm + '_4_60_mutation_codon_counts.df'))
    df_aas.to_pickle(join(General.Basepath, grpnm + '_4_60_mutation_aas_counts.df'))
    sms.to_pickle(join(General.Basepath, grpnm + '_4_60_nucleotide_count.df'))


