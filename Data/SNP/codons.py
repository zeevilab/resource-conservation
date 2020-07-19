'''
Genetic code table 11 - Bacteria, Archaea and Plant Plastid
'''
from pandas import DataFrame
from collections import defaultdict
from pandas.core.reshape.concat import concat

def get_codon_table(transl_table=11):
    if transl_table==11:
        codons = DataFrame(dict(
                AAs     = list('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
                Starts  = list('---M------**--*----M------------MMMM---------------M------------'),
                Base1   = list('TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'),
                Base2   = list('TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'),
                Base3   = list('TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG')
                ))
    elif transl_table==1:
        codons = DataFrame(dict(
                AAs     = list('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
                Starts  = list('---M------**--*----M---------------M----------------------------'),
                Base1   = list('TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'),
                Base2   = list('TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'),
                Base3   = list('TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG')
                ))
    else:
        raise ValueError('Table {} not coded'.format(transl_table))
    codons['Codon'] = codons[['Base1','Base2','Base3']].sum(1)
    codons = codons.set_index('Codon')
    starts = codons[codons.Starts == 'M'].dropna().index
    aas = codons['AAs']
    
    def determine_synonymity(codon, col = 'AAs'):
        ret = defaultdict(dict)
        for i, base in enumerate(codon.name):
            s = 0
            for opt in ['T','C','A','G']:
                if opt == base: continue
                newcodon = list(codon.name)
                newcodon[i] = opt
                newcodon = ''.join(newcodon)
                if codons.loc[newcodon][col] == codon[col]:
                    s += 1
            ret[codon.name][('S', i)] = s
            ret[codon.name][('NS', i)] = 3-s
        return DataFrame(ret).T
    
    synonymity = concat([determine_synonymity(codon) for nm, codon in codons.iterrows()])
    startsyn = concat([determine_synonymity(codon, 'Starts') for nm, codon in codons.iterrows() \
                       if codon.name in starts])
    return aas, startsyn, synonymity 