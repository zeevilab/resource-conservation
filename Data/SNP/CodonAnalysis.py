from pandas.io.pickle import read_pickle
from pandas.core.frame import DataFrame
from os.path import exists
from pandas.io.parsers import read_csv
from Data.SNP.codons import get_codon_table

def calc_codon_costs(out_f=None, force_rerun=False):
    if exists(out_f) and not force_rerun:
        return read_pickle(out_f)
    aas, _, _ = get_codon_table()
    ret = {}
    for cod_1, aa1 in aas.to_dict().items():
        for cod_2, aa2 in aas.to_dict().items():
            ret[(cod_1, cod_2)] = {'aa_s':aa1, 'aa_e':aa2}
    aa_props = read_csv('./aa_NCHP.csv', index_col=0)
    ret = DataFrame(ret).T
    ret.index.names = ['Codon_s','Codon_e']
    codon_props = ret.join(aa_props, on='aa_s').join(aa_props, on='aa_e', lsuffix='_s', rsuffix='_e')
    for v in ['C','N','hyd','PR']:
        codon_props[v+'_d'] = codon_props[v+'_e']-codon_props[v+'_s']
        codon_props[v+'_abs_d'] = abs(codon_props[v+'_d'])
    if out_f is not None:
        codon_props.to_pickle(out_f)
    return codon_props
