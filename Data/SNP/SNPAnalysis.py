import pysam
from operator import attrgetter
from itertools import groupby, product
from collections import defaultdict
from os.path import basename, join, exists
import numpy as np
from pandas.core.frame import DataFrame
from pandas.core.reshape.concat import concat
from pandas.io.pickle import read_pickle
from Bio import SeqIO
from scipy.stats import multinomial
from pandas.core.indexes.multi import MultiIndex
import logging
from collections import Counter
from pandas.core.series import Series
from lib import Utils
from Data.SNP.codons import get_codon_table
from Data.SNP.CodonAnalysis import calc_codon_costs
from Data.config import SNP
log_ = logging.getLogger(__name__)
del logging

def _getgenevariants(gene_name, var_group, only_snps, qual_thresh):
    lret = defaultdict(dict)
    for var in var_group:
        # Make sure variant has one letter in each allele.
        non_snp = False
        for allele in var.alleles:
            if len(allele) > 1:
                non_snp = True
                break
        if non_snp and only_snps:
            continue
        # Enforce quality threshold
        if qual_thresh is not None and var.qual < qual_thresh:
            continue
        for sample in var.samples.keys():
            val = var.samples[sample]['AD']
            for i, allele in enumerate(var.alleles):
                lret[(gene_name, var.pos, allele)][basename(sample).split('.')[0]] = val[i] \
                    if val[0] is not None else np.nan
    return DataFrame(lret)

def getvariants(vcf_f, outfol, only_snps, qual_thresh, min_samples, min_variants):
    varf = pysam.VariantFile(vcf_f)  # @UndefinedVariable
    ret = []
    for nm, grp in groupby(varf, attrgetter('contig')):
#         print(nm)
        lgrp = list(grp)
        # for each one, if total coverage is above something and any single variant depth is above something, report.
        genevars = _getgenevariants(nm, lgrp, only_snps, qual_thresh)
        # Min number of samples
        genevars = genevars[genevars.sum(1) > 0]
        if genevars.shape[0] < min_samples:
            continue
        # Min number of variants
        if genevars.groupby(level = 1, axis = 1).first().shape[1] < min_variants:
            continue
        ret.append(genevars.T)
    if len(ret) > 0:
        concat(ret, sort=False).to_pickle(join(outfol, '{}.df')\
                                          .format(basename(vcf_f.replace('.vcf',''))))
    
def _subsample(genepos, min_pos_reads):
    sampleto = genepos.sum()[genepos.sum() > min_pos_reads].min() # Note > and not >= for methods
    if np.isnan(sampleto):
        return genepos
    def samp(x):
        if x.sum() < sampleto:
            return x
        return multinomial.rvs(sampleto, (x/x.sum()).values)
    return genepos.fillna(0).apply(samp)
        
def _getgeneseqs(genes_df_f, db_fasta, gene_names, cachedir):
    cache_f = join(cachedir, basename(genes_df_f).replace('.df','.genes.dat'))
    if exists(cache_f):
        return Utils.Load(cache_f)
    ret = {}
    for rec in SeqIO.parse(db_fasta, 'fasta'):
        if rec.id in gene_names:
            ret[rec.id] = str(rec.seq)
            if len(ret) == len(gene_names):
                break
    Utils.Write(cache_f, ret)
    return ret

def _prepgene(nm, genedf, min_pos_reads, min_perc_poss, min_total_var_support, min_maf,
              min_samples, min_variants):
    # Calculate min number of positions for gene to be considered
    num_poss = genedf.groupby(level=[0,1]).first().shape[0]
    min_poss = max(num_poss*min_perc_poss/100., min_variants)
    # Take only samples that have the minimal number of reads per position for the minimal positions
    genedf = genedf.loc[:,(genedf.groupby(level=[0,1]).sum() >= min_pos_reads).sum() >= min_poss]
    # Discard if gene has less than minimal number of samples
    if genedf.shape[1] < min_samples:
        return
    # Remove positions that are not represented sufficiently in any sample
    ixs = (genedf.groupby(level = [0,1]).apply(lambda x:x.sum()>min_pos_reads).sum(1) == 0)\
                .replace(False, np.nan).dropna().index
    genedf = genedf[~genedf.index.get_level_values(1).isin(ixs)]
    # Make sure each variant has required support across samples and minimal maf to call SNP
    genedf = genedf[genedf.sum(1) > min_total_var_support]
    genedf = genedf[genedf.sum(1).groupby(level = 1)\
                            .apply(lambda x:(x.truediv(x.sum()) if x.sum() > 0 else 0)) > min_maf]
    genedf = genedf[genedf.sum(1).groupby(level = 1)\
                            .apply(lambda x:(x.truediv(x.sum()) if x.sum() > 0 else 0)) < 1-min_maf]
    remngposs = (genedf.sum(1).groupby(level = 1).count() > 1).replace(False, np.nan).dropna().index
    genedf = genedf.loc[genedf.index.get_level_values(1).isin(remngposs)]
    if genedf.groupby(level=1).first().shape[0] < min_variants:
        return
    #Subsampling
    genedf = genedf.groupby(level = 1).apply(_subsample, min_pos_reads)
    return genedf

def _get_consensus(genedf, geneseq):
    listseq = list(geneseq)
    for pos, snp in  genedf.sum(1).groupby(level = 1).apply(lambda x: x.idxmax()[2]).to_dict().items():
        listseq[pos-1] = snp
    return ''.join(listseq)
    
def _get_num_sites(codons, synonymity, starts, ignore_starts=False):
    if not ignore_starts and codons[0] in starts.index:
        sites = synonymity.loc[codons[1:]].sum().groupby(level=0).sum()
        sites = sites.add(starts.loc[codons[0]].groupby(level=0).sum())
        return sites
    return synonymity.loc[codons].sum().groupby(level=0).sum()

def calc_pnps(genedf, geneseq):
    aas, starts, synonymity = get_codon_table()
    consensus = _get_consensus(genedf, geneseq)
    if consensus is None:
        return None, None
    cons_codons = [consensus[i:i+3] for i in range(0,len(consensus),3)]
    num_sites = _get_num_sites(cons_codons, synonymity, starts)
    normdf = genedf.groupby(level=1).apply(lambda x:x.truediv(x.sum(0)))
    def determine_syn(x):
        oldcodon = cons_codons[(x.name[1]-1)//3]
        basepos = (x.name[1]-1)%3
        newbase = x.name[2]
        if oldcodon[basepos] == x.name[2]:
            return 0 # No mutation
        newcodon = oldcodon[:basepos] + newbase + oldcodon[basepos+1:]
        if (x.name[1]-1)//3 == 0 and oldcodon in starts.index:
            if newcodon in starts.index:
                return 'S'
            return 'NS'
        try:
            if aas[oldcodon] == aas[newcodon]:
                return 'S'
            return 'NS'
        except KeyError:
            return 0 # Can't determine synonymity
    normdf['Syn'] = normdf.apply(determine_syn, axis=1)
    sites = normdf.groupby('Syn').sum().loc[['NS','S']].fillna(0)
    dvalues = sites.truediv(num_sites, axis = 0)
    pnps = dvalues.loc['NS'].truediv(dvalues.loc['S'])
    sites['GeneSites'] = num_sites
    return sites, pnps

def get_groups(g1var1, g1val1, g1var2, g1val2, g2var1, g2val1, g2var2, g2val2, 
               g2_mirror_g1, g1_symmetrical):
    # Here we calculate the codon costs or just load them from file. 
    aa_props = calc_codon_costs(SNP.CodonProps)
    # Applying all the conditionals:
    aa_props = aa_props.dropna(subset = [g1var1] + ([g2var1] if g2var1 is not None else []) + \
                               ([g1var2] if g1var2 is not None else []) + \
                               ([g2var2] if g2var2 is not None else []))
    aa_g1 = aa_props[(aa_props[g1var1].apply(lambda x: eval(str(x) + g1val1))) &\
                       ((aa_props[g1var2].apply(lambda x: eval(str(x) + g1val2))) \
                        if g1var2 is not None else True)]
    if g1_symmetrical:
        aa_g1_1 = aa_props[(aa_props[g1var1].apply(lambda x: eval(str(x) + g1val2))) &\
                       ((aa_props[g1var2].apply(lambda x: eval(str(x) + g1val1))) \
                        if g1var2 is not None else True)]
        aa_g1 = concat([aa_g1, aa_g1_1], sort=False)
        aa_g1 = aa_g1.reset_index().drop_duplicates().set_index(['Codon_s','Codon_e'])
    if g2_mirror_g1:
        aag1ix = aa_g1.index.to_list()
        ixs = [ix for ix in aa_props.index.to_list() if ix not in aag1ix]
        aa_g2 = aa_props.loc[ixs]
    elif g2var1 is not None:
        aa_g2 = aa_props[(aa_props[g2var1].apply(lambda x: eval(str(x) + g2val1))) &\
                       ((aa_props[g2var2].apply(lambda x: eval(str(x) + g2val2))) \
                        if g2var2 is not None else True)]
    else:
        raise ValueError('Need either values or mirror for g2')
    # Remove synonymous mutations
    aa_g1 = aa_g1[aa_g1['aa_s'] != aa_g1['aa_e']]
    aa_g2 = aa_g2[aa_g2['aa_s'] != aa_g2['aa_e']]
    # Calculate the number of sites per codon
    ret = {}
    for p1 in ['A','C','G','T']:
        for p2 in ['A','C','G','T']:
            for p3 in ['A','C','G','T']:
                codon_s = p1+p2+p3
                ret[codon_s] = {('G1',0):0, ('G1',1):0, ('G1',2):0, ('G2',0):0, ('G2',1):0, ('G2',2):0}
                for mut_pos in [0,1,2]:
                    for mut_lett in ['A','C','G','T']:
                        codon_e = codon_s[:mut_pos] + mut_lett + codon_s[mut_pos+1:]
                        if codon_s == codon_e:
                            continue
                        if (codon_s, codon_e) in aa_g1.index:
                            ret[codon_s][('G1', mut_pos)] += 1
                        if (codon_s, codon_e) in aa_g2.index:
                            ret[codon_s][('G2', mut_pos)] += 1
    return set(aa_g1.index.to_list()), set(aa_g2.index.to_list()), DataFrame(ret).T

def create_all_pnpn_groups():
    # Creates all the groups we need once
    ret = {}
    # We iterate over all of the conditions we defined
    for g1var1, g1val1, g1var2, g1val2, g2var1, g2val1, g2var2, g2val2, g2mirrorg1, g1symm in \
        [['hyd_s', '<-2', 'hyd_e', '>2', None, None, None, None, True, True],
         ['hyd_s', '<-1', 'hyd_e', '>1', None, None, None, None, True, True],
         ['hyd_abs_d', '>1', None, None, 'hyd_abs_d', '<1', None, None, False, False],
         ['hyd_abs_d', '>2', None, None, 'hyd_abs_d', '<1', None, None, False, False],
         ['hyd_abs_d', '>2', None, None, 'hyd_abs_d', '<2', None, None, False, False],
         ['hyd_abs_d', '>3', None, None, 'hyd_abs_d', '<1', None, None, False, False],
         ['hyd_abs_d', '>3', None, None, 'hyd_abs_d', '<2', None, None, False, False],
         ['hyd_abs_d', '>3', None, None, 'hyd_abs_d', '<3', None, None, False, False],
         ['PR_abs_d', '>1', None, None, 'PR_abs_d', '<1', None, None, False, False],
         ['PR_abs_d', '>2', None, None, 'PR_abs_d', '<1', None, None, False, False],
         ['PR_abs_d', '>2', None, None, 'PR_abs_d', '<2', None, None, False, False],
         ['PR_abs_d', '>3', None, None, 'PR_abs_d', '<1', None, None, False, False],
         ['PR_abs_d', '>3', None, None, 'PR_abs_d', '<2', None, None, False, False],
         ['PR_abs_d', '>3', None, None, 'PR_abs_d', '<3', None, None, False, False]]:
        # Construct the name of each condition for the dataframe
        nm = 'G1:{}{}'.format(g1var1,g1val1)
        nm += '&{}{}'.format(g1var2,g1val2) if g1var2 is not None else ''
        nm += '|{}{}&{}{}'.format(g1var1,g1val2,g1var2,g1val1) if g1symm else ''
        nm += '__G2:{}{}'.format(g2var1,g2val1) if g2var1 is not None else ''
        nm += '&{}{}'.format(g2var2,g2val2) if g2var2 is not None else ''
        nm += '__G2:mirror' if g2mirrorg1 else ''
        g1, g2, nononymity = get_groups(g1var1, g1val1, g1var2, g1val2, g2var1, g2val1, g2var2, g2val2, 
                                       g2mirrorg1, g1symm)
        ret[nm] = [g1,g2,nononymity]
    return ret

def calc_pnpn_one(cons_codons, normdf, g1, g2, nononymity):
    # Get the number of sites for each group for this gene
    num_sites = _get_num_sites(cons_codons, nononymity, [], True)
    # This internal method determines whether a SNP belongs to group 1 (G1), group 2 (G2), 
    # neither group (G0) or has no mutation (0). Only G1 and G2 are subsequently used.
    def determine_nyn(x):
        oldcodon = cons_codons[(x.name[1]-1)//3]
        basepos = (x.name[1]-1)%3
        newbase = x.name[2]
        if oldcodon[basepos] == x.name[2]:
            return 0 # No mutation
        newcodon = oldcodon[:basepos] + newbase + oldcodon[basepos+1:]
        if (oldcodon, newcodon) in g1:
            return 'G1'
        elif (oldcodon, newcodon) in g2:
            return 'G2'
        else:
            return 'G0'
    # Here we call the above method for each sample
    normdf['Nyn'] = normdf.apply(determine_nyn, axis=1)
    # Sum up for every sample and add the number of sites per this gene.
    sites = normdf.groupby('Nyn').sum().reindex(['G1','G2']).fillna(0)
    sites['GeneSites'] = num_sites
    return sites
    
def calc_pnpn_all(genedf, geneseq, pnpn_groups):
    # LIAT:
    # Here we calculate pN(group1)/pN(group2) for all different groups that we define below.
    # We start by taking the consensus sequence:
    consensus = _get_consensus(genedf, geneseq)
    if consensus is None:
        return None
    # Splitting to codons
    cons_codons = [consensus[i:i+3] for i in range(0,len(consensus),3)]
    # Separately we normalize by the numebr of reads in each poisition to get a relative abundance
    # of each SNP
    normdf = genedf.groupby(level=1).apply(lambda x:x.truediv(x.sum(0)))
    ret = []
    for nm, (g1, g2, nononymity) in pnpn_groups.items():
        # Call a calculation for each condition
        res = calc_pnpn_one(cons_codons, normdf, g1, g2, nononymity)
        # Organize with the constructed name as index
        res.index = MultiIndex.from_tuples([(nm, i) for i in res.index])
        ret.append(res)
    return concat(ret, sort=False)

def _get_ffdeg_sites(cons_codons, synonymity):
    ff_syn_codons = set(synonymity[synonymity[('S',2)] == 3].index)
    poss = [(i+1)*3 for i,cod in enumerate(cons_codons) if cod in ff_syn_codons]
    return poss

def calc_ffdeg_pi_within(genedf, geneseq):
    '''
    pi(S,G) = 1/|Gffdeg| \sum_(i=1)^|Gffdeg| \sum_(B_1\in{ACTG}) \sum_(B_2\in{ACTG}-B_1) 
                (X_(i,B_1) * X_(i,B_2)) / (C_i * (C_i - 1))
    where S is the sample, G is the gene of interest, |Gffdeg| is the number of fourfold degenerate 
    sites in the gene, X_(i,Bj) the number of nucleotide Bj seen at fourfold degenerate position i 
    and c_i the coverage at position i
    '''
    _, _, synonymity = get_codon_table()
    consensus = _get_consensus(genedf, geneseq)
    if consensus is None:
        return None
    cons_codons = [consensus[i:i+3] for i in range(0,len(consensus),3)]
    ffdeg_poss = _get_ffdeg_sites(cons_codons, synonymity)
    genedf_filt = genedf[genedf.index.get_level_values(1).isin(ffdeg_poss)]
    def calconepos(pos):
        tosum = []
        cov = pos.sum()
        for (nm1, XiB1), (nm2, XiB2) in product(pos.iterrows(), pos.iterrows()):
            if nm1[2] == nm2[2]: continue
            tosum.append(XiB1.multiply(XiB2))
        return (sum(tosum)).truediv(cov.multiply(cov-1))
    return genedf_filt.groupby(level = 1).apply(calconepos).sum(0).truediv(len(ffdeg_poss)), len(ffdeg_poss)

def calc_percodon_mut(genedf, geneseq):
    consensus = _get_consensus(genedf, geneseq)
    if consensus is None:
        return None, None
    cons_codons = [consensus[i:i+3] for i in range(0,len(consensus),3)]
    num_sites = Counter(cons_codons)
    normdf = genedf.groupby(level=1).apply(lambda x:x.truediv(x.sum(0)))
    def count_muts(x):
        oldcodon = cons_codons[(x.name[1]-1)//3]
        basepos = (x.name[1]-1)%3
        newbase = x.name[2]
        if oldcodon[basepos] == x.name[2]:
            return (oldcodon, oldcodon) # No mutation
        newcodon = oldcodon[:basepos] + newbase + oldcodon[basepos+1:]
        return (oldcodon, newcodon)
    normdf['MutDir'] = normdf.apply(count_muts, axis=1)
    sites = normdf.groupby('MutDir').sum()
    sites = sites.loc[[i for i in sites.index if i[0] != i[1]]]  
    return [sites, num_sites]
                                         
def analyze_one(nm, genedf, min_pos_reads, min_perc_poss, min_total_var_support, min_maf,
                  min_samples, min_variants, gene_seqs, pnpn_groups):
    # This method prepares a separate dataframe for each gene (subsampling etc.) and then calculates
    # all of the required metrics.  
    genedf = _prepgene(nm, genedf, min_pos_reads, min_perc_poss, min_total_var_support, min_maf,
                       min_samples, min_variants)
    if genedf is None:
        return [None] * 5
    sites, pnps = calc_pnps(genedf, gene_seqs[nm])
    percmut = calc_percodon_mut(genedf, gene_seqs[nm])
    pnpn = calc_pnpn_all(genedf, gene_seqs[nm], pnpn_groups)
    pnpn.index = MultiIndex.from_tuples([(nm,ix[0],ix[1]) for ix in pnpn.index])
    for df in ([sites, percmut[0]] if sites is not None else [percmut[0]]):
        df.index = MultiIndex.from_product([[nm], df.index])
    percmut[1] = Series(percmut[1]).to_frame(nm).T
    ffdeg_pw, ffdeg_poss = calc_ffdeg_pi_within(genedf, gene_seqs[nm])
    return sites, (pnps.to_frame(nm).T if pnps is not None else pnps), percmut, \
            (ffdeg_pw.to_frame(nm).T, ffdeg_poss), pnpn
    

def analyze_genes(genes_df_f, db_fasta, outdir, cachedir, min_pos_reads, min_perc_poss, 
                  min_total_var_support, min_maf, min_samples, min_variants):
    # This method potentially calculates all metrics (pN/pS, pi_within, etc.) for one gene file
    # out of approximately 2500, containing up to 10,000 separate genes (usually 500-1000)
    # It iterates all genes and calculates everything for each.
    # Then it concatenates and saves everything. The calling for each gene is to 'analyze_one'. 
    df = read_pickle(genes_df_f)
    pnpn_groups = create_all_pnpn_groups()
    log_.info('Analyzing file {} minpos {} minperc {}'.format(basename(genes_df_f), min_pos_reads, 
                                                              min_perc_poss))
    gene_seqs = _getgeneseqs(genes_df_f, db_fasta, df.index.get_level_values(0).unique(), cachedir)
    sitess, pnpss, percmuts, ffdeg_pws, pnpns = [], [], [], [], []
    ffdeg_poss = {}
    for nm, genedf in df.groupby(level=0):
        sites, pnps, percmut, ffdeg_pw, pnpn  = analyze_one(nm, genedf, 
                                                        min_pos_reads, min_perc_poss, 
                                                        min_total_var_support, min_maf, min_samples, 
                                                        min_variants, gene_seqs, pnpn_groups)
        if sites is not None:
            sitess.append(sites)
        if pnps is not None:
            pnpss.append(pnps)
        if percmut is not None:
            percmuts.append(percmut)
        if ffdeg_pw is not None:
            ffdeg_pws.append(ffdeg_pw[0])
            ffdeg_poss[nm] = ffdeg_pw[1]
        if pnpn is not None:
            pnpns.append(pnpn)
    baseout = basename(genes_df_f).replace('.df','')
    if len(sitess) > 0 and any([s is not None for s in sitess]):
        concat(sitess, sort=False).to_pickle(join(outdir, baseout + '.muts.df'))
    if len(pnpss) > 0 and any([p is not None for p in pnpss]):
        concat(pnpss, sort=False).to_pickle(join(outdir, baseout + '.pnps.df'))
    if len(percmuts) > 0:
        concat([p[0] for p in percmuts], sort=False).to_pickle(join(outdir, baseout + '.percmut.df'))
        concat([p[1] for p in percmuts], sort=False).to_pickle(join(outdir, baseout + '.codons.df'))
    if len(ffdeg_pws) > 0:
        concat(ffdeg_pws, sort=False).to_pickle(join(outdir, baseout + '.ffdeg_pi_wit.df'))
        Series(ffdeg_poss).to_pickle(join(outdir, baseout + '.ffdeg_poss.df'))
    if len(pnpns) > 0:
        concat(pnpns, sort=False).to_pickle(join(outdir, baseout + '.pnpn.df'))
