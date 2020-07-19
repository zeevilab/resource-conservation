from os.path import join, exists
from pandas.io.pickle import read_pickle
from pandas.core.reshape.concat import concat
import numpy as np
from _collections import defaultdict
from glob import glob
from pandas.core.frame import DataFrame
from builtins import dict
from lib.Utils import mkdirifnotexists, split3way
from lib import Utils
from Data.config import SNP, General

def do_one_group_pnps(nm, f_prefix, genegroup, analysisclass, mingenes=3):
    ret = []
    for prefix in set([g[:-4] for g in genegroup]):
        f_in = join(analysisclass.OutDir, prefix + '.muts.df')
        if not exists(f_in):
            continue
        try:
            ldf = read_pickle(f_in)
            ldf = ldf.loc[[g for g in genegroup if g in ldf.index]]
        except KeyError:
            continue
        if ldf.shape[0] > 0:
            ret.append(ldf)
    if len(ret) == 0:
        return
    outdf = concat(ret, sort = False)
    outpath = mkdirifnotexists(join(analysisclass.OutDirCollate, 'pnps')) 
    if outdf.groupby(level=0).first().shape[0] < mingenes:
        return
    outdf.to_pickle(join(outpath, f_prefix + '_' + nm + '.pnps.df'))
    
def do_one_group_pnpn(nm, f_prefix, genegroup, analysisclass, mingenes=3):
    ret = []
    for prefix in set([g[:-4] for g in genegroup]):
        f_in = join(analysisclass.OutDir, prefix + '.pnpn.df')
        if not exists(f_in):
            continue
        try:
            ldf = read_pickle(f_in)
            ldf = ldf.loc[[g for g in genegroup if g in ldf.index]] 
        except KeyError:
            continue
        if ldf.shape[0] > 0:
            ret.append(ldf)
    if len(ret) == 0:
        return
    outdf = concat(ret, sort = False)
    outpath = mkdirifnotexists(join(analysisclass.OutDirCollate, 'pnpn'))
    if outdf.groupby(level=0).first().shape[0] < mingenes:
        return
    outdf.to_pickle(join(outpath, f_prefix + '_' + nm + '.pnpn.df'))
    
def do_one_group_ffdeg_piwit(nm, f_prefix, genegroup, analysisclass, mingenes=3):
    ret = []
    retlens = []
    for prefix in set([g[:-4] for g in genegroup]):
        f_in = join(analysisclass.OutDir, prefix + '.ffdeg_pi_wit.df')
        if not exists(f_in):
            continue
        try:
            ldf = read_pickle(f_in)
            ldf = ldf.loc[[g for g in genegroup if g in ldf.index]]
            llensdf = read_pickle(f_in.replace('_pi_wit','_poss'))
            llensdf = llensdf.loc[[g for g in genegroup if g in llensdf.index]]
        except KeyError:
            continue
        if ldf.shape[0] > 0:
            ret.append(ldf)
            retlens.append(llensdf)
    if len(ret) == 0:
        return
    outdf = concat(ret, sort = False).dropna(how='all').dropna(how='all', axis = 1)
    if outdf.shape[0] < mingenes:
        return
    outdf_lens = concat(retlens)
    outdf_lens.name = 'Length'
    ret = {}
    for col in outdf:
        coldf = outdf[[col]].multiply(outdf_lens,axis=0).join(outdf_lens).dropna().sum()
        ret[col] = {'pi':coldf[col] / coldf['Length'], 'length':coldf['Length'], 
                    'num_genes':len(outdf[col].dropna())}
    outpath = mkdirifnotexists(join(analysisclass.OutDirCollate, 'ffdeg'))
    DataFrame(ret).to_pickle(join(outpath, f_prefix + '_' + nm + '.ffdeg_pi_wit.df'))

def filter_db(dbdct, analysisclass, mingenes):
    indir = analysisclass.OutDir
    pnps_piwig_f = join(indir, 'PNPSPiWiGenes.dat')
    genes_in_use = []
    if exists(pnps_piwig_f):
        pnpsgenes = Utils.Load(pnps_piwig_f)
    else:
        for fname in glob(join(indir, '*pnps.df')) + glob(join(indir, '*pi_wit.df')):
            print(fname)
            genes_in_use.extend(list(read_pickle(fname).index.get_level_values(0).unique()))
        pnpsgenes = list(set(genes_in_use))
        Utils.Write(pnps_piwig_f, genes_in_use)
    ret = {}
    pnpsgenes = set(pnpsgenes)
    for k,v in dbdct.items():
        if len(set(v).intersection(pnpsgenes)) >= mingenes:
            ret[k] = v
    return(ret)
    
def do_collate(f_prefixes, minpos, minperc, mingenes, minsamples, minsamples_gene):
    ret = defaultdict(dict)
    ps = defaultdict(dict)
    pn = defaultdict(dict)
    for fname in glob(f_prefixes + '*.pnps.df'):
        grpname = split3way(fname)[1].replace('.pnps','')
        df = read_pickle(fname)
        keepinds = df.index.get_level_values(0).isin(\
                ((df.groupby(level=0).count() > 1).sum(1) >= minsamples_gene)\
                .replace(False, np.nan).dropna().index)
        df = df.loc[keepinds]
        df = df.loc[:,(df.groupby(level=0).count() > 1).sum(0) >= mingenes]
        if df.shape[1] <= minsamples:
            continue
        gs = df[['GeneSites']]
        df = df.drop('GeneSites', axis = 1)
        for col in df.columns:
            coldf = df[[col]].join(gs).dropna()
            coldf = coldf.groupby(level=1).sum()
            coldf = coldf[col].truediv(coldf['GeneSites'])
            ret[grpname][col] = coldf.NS/coldf.S
            pn[grpname][col] = coldf.NS
            ps[grpname][col] = coldf.S
        print(grpname)
    df = DataFrame(ret)
    ps = DataFrame(ps)
    pn = DataFrame(pn)
    df.to_csv(join(General.Basepath, '{}_{}_{}_{}_{}_{}.csv'\
                                  .format(f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))
    pn.to_csv(join(General.Basepath, '{}_{}_{}_{}_{}_{}.pn.csv'\
                                  .format(f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))
    ps.to_csv(join(General.Basepath.OutDirCollate, '{}_{}_{}_{}_{}_{}.ps.csv'\
                                  .format(f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))

def _collate_pnpn_inner(fname, mingenes, minsamples, minsamples_gene, outdir):
    grpname = split3way(fname)[1].replace('.pnpn','').split(':')[1]
    ret = defaultdict(dict)
    ret_g1 = defaultdict(dict)
    ret_g2 = defaultdict(dict)
    df = read_pickle(fname)
    for nm, ldf in df.groupby(level=1):
        keepinds = ldf.index.get_level_values(0).isin(\
                ((ldf.groupby(level=0).count() > 1).sum(1) >= minsamples_gene)\
                .replace(False, np.nan).dropna().index)
        ldf = ldf.loc[keepinds]
        ldf = ldf.loc[:,(ldf.groupby(level=0).count() > 1).sum(0) >= mingenes]
        if ldf.shape[1] <= minsamples:
            continue
        gs = ldf[['GeneSites']]
        ldf = ldf.drop('GeneSites', axis = 1)
        for col in ldf.columns:
            coldf = ldf[[col]].join(gs).dropna()
            coldf = coldf.groupby(level=2).sum()
            coldf = coldf[col].truediv(coldf['GeneSites'])
            ret_g1[(nm,grpname)][col] = coldf['G1']
            ret_g2[(nm,grpname)][col] = coldf['G2']
            ret[(nm,grpname)][col] = (coldf['G1']/coldf['G2']) if coldf['G2'] !=0 else np.nan
    outdf = DataFrame(ret)
    outdf_g1 = DataFrame(ret_g1)
    outdf_g2 = DataFrame(ret_g2)
    if outdf.shape != (0,0):
        outdf.to_pickle(join(outdir, grpname + '.tmp.df'))
    if outdf_g1.shape != (0,0):
        outdf_g1.to_pickle(join(outdir, grpname + '.tmp.g1.df'))
    if outdf_g2.shape != (0,0):
        outdf_g2.to_pickle(join(outdir, grpname + '.tmp.g2.df'))
    
def do_collate_pnpn(f_prefixes, minpos, minperc, mingenes, minsamples, minsamples_gene):
    tmpdir = mkdirifnotexists(join(SNP.OM_RGC.OutDirCollate, 'pnpn', 'tmpfiles'))
    for fname in glob(f_prefixes + '*.pnpn.df'):
        _collate_pnpn_inner(fname, mingenes, minsamples, minsamples_gene, tmpdir)
    ret = []
    ret_g1 = []
    ret_g2 = []
    for fname in glob(join(tmpdir, '*.tmp.df')):
        ret.append(read_pickle(fname).T)
    for fname in glob(join(tmpdir, '*.tmp.g1.df')):
        ret_g1.append(read_pickle(fname).T)
    for fname in glob(join(tmpdir, '*.tmp.g2.df')):
        ret_g2.append(read_pickle(fname).T)
    outdir = join(SNP.OM_RGC.OutDirCollate, 'pnpn')
    with open(join(outdir, 'pNpNCases.txt'), 'w') as ftxt:
        ftxt.write('Conditions for pN groups in this analysis\n')
        ftxt.write('Always pN(G1)/pN(G2) so invert if G1 is more conservative\n\n')
        bigdf = concat(ret, sort=False)
        bigg1 = concat(ret_g1, sort=False)
        bigg2 = concat(ret_g2, sort=False)
        for j, col in enumerate(bigdf.index.get_level_values(0).unique()):
            ftxt.write('Case {}: {}\n'.format(j,col))
            bigdf.loc[col].T.to_csv(join(SNP.OM_RGC.OutDirCollate, 'pnpn', 
                                            'pNpN_Case_{}_{}_{}_{}_{}_{}_{}.csv'\
                                  .format(j, f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))
            bigg1.loc[col].T.to_csv(join(SNP.OM_RGC.OutDirCollate, 'pnpn', 
                                            'pNpN_Case_{}_{}_{}_{}_{}_{}_{}.g1.csv'\
                                  .format(j, f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))
            bigg2.loc[col].T.to_csv(join(SNP.OM_RGC.OutDirCollate, 'pnpn', 
                                            'pNpN_Case_{}_{}_{}_{}_{}_{}_{}.g2.csv'\
                                  .format(j, f_prefixes.split('/')[-1], minpos, minperc, mingenes, 
                                          minsamples, minsamples_gene)))
    
