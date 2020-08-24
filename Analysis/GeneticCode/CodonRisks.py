from os.path import join
import numpy as np
from pandas.core.frame import DataFrame
from pandas.io.parsers import read_csv
from collections import defaultdict
from Data.config import CodeAnalysis
from Analysis.GeneticCode.CodonTransitions import get_relevant_codonpairs
from lib.Utils import mkdirifnotexists
from lib import Utils

def determine_mutpos(x):
    x0, x1 = x.name
    mutpos =  [i for i in range(len(x0)) if x0[i] !=x1[i]]
    if len(mutpos) == 0:
        return -1
    return mutpos[0]

def determine_muttype(x, all_mutations):
    pos = int(x.mutation_pos)
    if pos == -1:
        return 'None'
    x0, x1 = x.name
    if not all_mutations:
        if (x0[pos] in 'AG' and x1[pos] in 'AG') or (x0[pos] in 'CT' and x1[pos] in 'CT'):
            return 'Transition'
        return 'Transversion'
    else:
        return ''.join([x0[pos],x1[pos]])
    
def filter_nonsyn(mut_counts, aas, only_fourfold_deg=False):
    mut_counts['aa_start'] = [aas[ix] for ix in mut_counts.index.get_level_values(0)]
    mut_counts['aa_end'] = [aas[ix] for ix in mut_counts.index.get_level_values(1)]
    if only_fourfold_deg:
        def isffdeg(x):
            xff = x[x.index.get_level_values(1).str.startswith(x.name[:2])]
            return (xff.aa_start == xff.aa_end).all()
        ff_deg = mut_counts.groupby(level=0).apply(isffdeg).replace(False,np.nan).dropna().index
        mut_counts = mut_counts.loc[ff_deg].groupby(level=0).apply(lambda x: x.loc[x.name]\
                                .loc[x.index.get_level_values(1).str.startswith(x.name[:2])])
    return mut_counts[(mut_counts['aa_start']==mut_counts['aa_end']) \
                            & (mut_counts['aa_start']!='*')]\
                            .drop(['aa_start','aa_end'], axis=1)

def getmuts(df, all_mutations):
    df['mutation_pos'] = df.apply(determine_mutpos, axis=1)
    df['mutation_type'] = df.apply(determine_muttype, args=[all_mutations], axis=1)
    return df

def get_mut_costs(aas):
    codon_pairs = get_relevant_codonpairs()
    aa_muts = DataFrame({(i,j):{'aa_start':aas[i], 'aa_end':aas[j]} for i,j in codon_pairs}).T
    
    aa_props = read_csv('./resource/aa_s.csv').set_index('Code')
    aa_muts = aa_muts.join(aa_props.rename(columns={c:c+'_s' for c in aa_props.columns}), on='aa_start')\
                .join(aa_props.rename(columns={c:c+'_e' for c in aa_props.columns}), on='aa_end')
    for col in aa_muts.columns:
        if col.endswith('_s'):
            aa_muts[col.replace('_s','_d')] = aa_muts[col.replace('_s','_e')].subtract(aa_muts[col])
    aa_muts = aa_muts[[c for c in aa_muts.columns if c.endswith('_d')]]
    aa_muts.index.names = ['aa_start','aa_end']
    aa_muts.Hyd_d = aa_muts.Hyd_d.abs()
    aa_muts.PR_d = aa_muts.PR_d.abs()
    return aa_muts

def scramble_codons(aas):
    aa_starts = list(set([i[:2] for i in aas.index]))
    aa_trans = dict(zip(aa_starts, np.random.permutation(aa_starts)))
    # Preserving stop codons close in permutation
    stop_sq_1 = aa_trans['TA']
    stop_sq_2 = aa_trans['TG']
    def get_transition(x):
        return 'A' if x=='G' else 'G' if x=='A' else 'C' if x=='T' else 'T'
    if np.random.rand() > 0.5:
        desired_stop_sq_2 = ''.join([stop_sq_1[0], get_transition(stop_sq_1[1])])
    else:
        desired_stop_sq_2 = ''.join([get_transition(stop_sq_1[0]), stop_sq_1[1]])
    newcode = {}
    for k, v in aa_trans.items():
        if v==stop_sq_2: newcode[k] = desired_stop_sq_2
        elif v==desired_stop_sq_2: newcode[k] = stop_sq_2
        else: newcode[k] = v
    aa_trans = newcode
    aas_new = aas.copy()
    aas_new.index = [''.join([aa_trans[i[:2]],i[2]]) for i in aas.index]
    return aas_new

def codon_risk(df, aas, prefix, all_mutations=True, external_counts=None, 
               external_titv=None, subdir=None):
    stops = ['TAA','TAG','TGA']
    if external_counts is None:
        df = df.sum(1).to_frame('all_muts')
        dfnostop = df.loc[[(i,j) for i,j in df.index if i not in stops and j not in stops]]
        codonabuns = df.groupby(level=0).sum()
    else:
        codonabuns = external_counts.to_frame('all_muts')
    codonabuns['AAs'] = [aas[i] for i in codonabuns.index]
    meanaas = codonabuns.groupby('AAs').median()
    codonabuns = codonabuns.join(meanaas, on='AAs', lsuffix='_1').drop(['all_muts_1','AAs'], axis=1)
    codonabuns = codonabuns.truediv(codonabuns.sum()).rename(columns={'all_muts':'codon_abuns'})
    ret = defaultdict(list)
    for it in range(10001):
        if it == 0:
            newaas = aas
            newcodonabuns = codonabuns
        else:
            # To maintain the abundances of the codons coding for the same amino acids
            codonaas = codonabuns.copy()
            codonaas['AAs'] = [aas[i] for i in codonabuns.index]
            codonaas = codonaas.reset_index().groupby('AAs').apply(lambda x:x.reset_index())\
                        .drop(['level_0','AAs'], axis=1)
            newaas = scramble_codons(aas)
            codon_shuf = newaas.reset_index().groupby('AAs').apply(lambda x:x.reset_index())\
                        .drop(['level_0','AAs'], axis=1)
            newcodonabuns = codonaas.join(codon_shuf, lsuffix='_1').set_index('index')['codon_abuns']
        # Estimate the abundance of mutations using fourfold degenerate synonymous mutations
        if external_titv is None:
            mutabuns = getmuts(filter_nonsyn(dfnostop, aas, True), all_mutations)\
                        .drop('mutation_pos', axis=1).groupby('mutation_type').mean().drop('None')
            mutabuns = mutabuns.truediv(mutabuns.sum()).rename(columns={'all_muts':'mut_abuns'})
        else:
            mutabuns = DataFrame({'mut_abuns':{'Transition':external_titv[0],
                                               'Transversion':external_titv[1]}})
            mutabuns.index.name = 'mutation_type'
        mut_costs = get_mut_costs(newaas)
        newstops = newaas[newaas=='*'].index
        allabuns = getmuts(mut_costs, all_mutations).join(mutabuns, on='mutation_type').reset_index()\
                    .join(newcodonabuns, on='aa_start').set_index(['aa_start','aa_end'])\
                    .drop(['mutation_pos','mutation_type'], axis=1)
        allabuns = allabuns.loc[[(i,j) for i,j in allabuns.index \
                                 if i not in newstops and j not in newstops]]
        def applyfunc(x, col):
            return x[col]*x.mut_abuns*x.codon_abuns
        ret['hyd_risk'].append(allabuns.apply(lambda x:applyfunc(x, col='Hyd_d'), axis=1).sum())
        ret['PR_risk'].append(allabuns.apply(lambda x:applyfunc(x, col='PR_d'), axis=1).sum())
        ret['n+_risk'].append(allabuns[allabuns.N_d>0].apply(lambda x:applyfunc(x, col='N_d'), axis=1).sum())
        ret['c+_risk'].append(allabuns[allabuns.C_d>0].apply(lambda x:applyfunc(x, col='C_d'), axis=1).sum())
        ret['o+_risk'].append(allabuns[allabuns.O_d>0].apply(lambda x:applyfunc(x, col='O_d'), axis=1).sum())
        ret['code'].append(newaas)
    if subdir is not None: outfol = mkdirifnotexists(join(CodeAnalysis.CodonsDir, subdir))
    outfname = 'Codon_risk_{}_{}_pstop_medaa.dat'.format('allmut' if all_mutations else 'TiTv',
                                                   prefix)
    print(outfname)
    Utils.Write(join(outfol, outfname), ret)
    for k in ret:
        if k != 'code':
            print('{}: {}'.format(k, sum(ret[k]<ret[k][0])/10000.))
    
