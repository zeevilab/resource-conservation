import numpy as np
from itertools import product
import logging
import gem_qa
log_ = logging.getLogger(__name__)
from time import sleep
max_1 = lambda a: a if a>0 else 1

akmer_dict = {''.join(lst): {ord(bp):''.join(lst[1:] + (bp,)) for bp in ['A', 'G', 'C', 'T', 'N']}  for lst in product(['A', 'G', 'C', 'T', 'N'], repeat = 4)}
fmerdict = {''.join(tmer):i for i, tmer in enumerate(product(['A', 'C', 'G', 'T', 'N'], repeat = 4))}
NNNN = fmerdict['NNNN']
fmer_ar = [None] * 3125
for tmer in product(['A', 'C', 'G', 'T', 'N'], repeat = 4):
    for b_i, b in enumerate(['A', 'G', 'C', 'T', 'N']):
        fmer_ar[fmerdict[''.join(tmer)] * 5 + b_i] = fmerdict[''.join(tmer[1:] + (b, ))]

conv_arr = [None] * 30
x = 0
while x < 30:
    conv_arr[x] = 100
    x += 1
A_C = ord('A')
conv_arr[0] = 0 #'A'
conv_arr[6] = 1 #'G'
conv_arr[2] = 2 #'C'
conv_arr[19] = 3 #'T'
conv_arr[13] = 4 #'N'
conv_arr[24] = 4 #'Y' => 'N'
conv_arr[17] = 4 #'R' => 'N'
conv_arr[22] = 4 #'W' => 'N'
conv_arr[18] = 4 #'S' => 'N'
conv_arr[10] = 4 #'K' => 'N'
conv_arr[12] = 4 #'M' => 'N'
conv_arr[3] = 4 #'D' => 'N'
conv_arr[21] = 4 #'V' => 'N'
conv_arr[7] = 4 #'H' => 'N'
conv_arr[1] = 4 #'B' => 'N'
conv_arr[23] = 4 #'X' => 'N'

def conv_ar(c):
    return conv_arr[c - A_C]

def validate_fmerdict(fmerdict):
    for k in fmerdict.keys():
        assert fmerdict[k] == fmer_ar[fmer_ar[fmer_ar[fmer_ar[NNNN * 5 + conv_ar(ord(k[0]))] * 5 + conv_ar(ord(k[1]))] * 5 + conv_ar(ord(k[2]))] * 5 + conv_ar(ord(k[3]))]
        for b in 'AGCTN':
            assert fmerdict[k] == fmer_ar[fmerdict[b + k[:-1]] * 5 + conv_ar(ord(k[-1]))]

THRESH = 5000

def map_qual_score(seq, quals, gigar, dest_id, pos, strand, lseq, mismatches, dc,
                    mdl0, jnd0, mdl0s5, jnd0s5, mdl0s2, jnd0s2, mdl0s25, jnd0s25,
                    mdl0s12, jnd0s12, mdl0s125, jnd0s125, refins, refins_lgsum, rdins, rdins_lgsum):
    if strand == 1:
        curref = dc._get_seq(dest_id)
    else:
        curref = dc._get_rev(dest_id)
    len_ref = dc.get_len(dest_id)
    
    if '>' in gigar:
        if strand == 1:
            ref_of = pos - 1
            end = pos - 1 + lseq
        else:
            ref_of = len_ref - pos + 1 - lseq 
            end = len_ref - pos + 1
         
        ga, ofst = gem_qa.get_opt_gigar_array(len_ref, curref, ref_of%len_ref, end%len_ref, strand, gigar, lseq, seq, mismatches)
        lga = len(ga) 
    else:
        ofst = 0
        ga = gem_qa.create_gigar_array(gigar)
        lga = len(ga)
    if strand == 1:
        ref_of = pos - 1 + ofst
        end = pos - 1 + lseq + ofst + gem_qa.sum_skips(ga)
    else:
        ref_of = len_ref - pos  + 1 - lseq - ofst - gem_qa.sum_skips(ga)
        end = len_ref - pos + 1 - ofst 
    ref_of = ref_of % len_ref
    end = end % len_ref
    
    debug = 0
    prevmer = NNNN
    res = 0
    r_i = 0
    d_i = 0
    g_i = 0
    stupid = 1
    
    while g_i < lga:
        gb = ga[g_i]
        if gb == 1:
            inslen = 1
            while g_i+1 < lga and ga[g_i+1] == 1:
                g_i+=1
                inslen += 1
            res += refins[inslen] - refins_lgsum
            d_i += inslen
        elif gb == -1:
            inslen = 1
            prevmer = fmer_ar[prevmer * 5 + conv_arr[ord(seq[r_i]) - A_C]]
            r_i += 1
            while g_i+1 < lga and ga[g_i+1] == -1:
                g_i += 1
                inslen += 1
                prevmer = fmer_ar[prevmer * 5 + conv_arr[ord(seq[r_i]) - A_C]]
                r_i += 1
            res += rdins[inslen] - rdins_lgsum
        else:
            convbp = conv_arr[ord(seq[r_i]) - A_C]
            convrefbp = conv_arr[ord(curref[ref_of+d_i]) - A_C]
            convrefbp = convbp if convrefbp == 4 else convbp
            cq = quals[r_i]
            stupid = 1
            while stupid == 1:
                stupid = 0
                curtot_sum = mdl0s5[r_i, prevmer, cq, convrefbp]
                if curtot_sum > THRESH:
                    v = mdl0[r_i, prevmer, cq, convrefbp, convbp]
                    break
                
                curtot_sum = mdl0s25[r_i, cq, convrefbp]
                if curtot_sum > THRESH:
                    v = mdl0s2[r_i, cq, convrefbp, convbp]
                    break
                
                curtot_sum = jnd0s5[r_i, prevmer, cq, convrefbp]
                if curtot_sum > THRESH:
                    v = jnd0[r_i, prevmer, cq, convrefbp, convbp]
                    break
                
                curtot_sum = jnd0s25[r_i, cq, convrefbp]
                if curtot_sum > THRESH:
                    v = jnd0s2[r_i, cq, convrefbp, convbp]
                    break
                
                curtot_sum = mdl0s125[cq, convrefbp]
                if curtot_sum > THRESH:
                    v = mdl0s12[cq, convrefbp, convbp]
                    break
                
                curtot_sum = jnd0s125[cq, convrefbp]
                if curtot_sum > THRESH:
                    v = jnd0s12[cq, convrefbp, convbp]
                    break
                
                return -100
            
            res += np.log10(max_1(v)) - np.log10(curtot_sum)
            d_i += 1
            if ref_of + d_i >= len_ref:
                ref_of = (ref_of + d_i)%len_ref - d_i
            r_i += 1
            prevmer = fmer_ar[prevmer * 5 + convbp]
        g_i += 1
    return 10**res