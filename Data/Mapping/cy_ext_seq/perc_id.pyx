cimport cython
from libcpp.string cimport string
from libc.math cimport pow, log10
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import numpy as np
cimport numpy as np
from libc.string cimport memset, strrchr, strncpy, strcmp
from itertools import product
cimport gem_qa 
cimport CSeqDict
import logging
log_ = logging.getLogger(__name__)
FTYPE = np.float64
ctypedef np.float64_t FTYPE_t
from time import sleep
cdef inline FTYPE_t max_1(FTYPE_t a): return a if a >0 else 1

cdef dict akmer_dict = {''.join(lst): {ord(bp):''.join(lst[1:] + (bp,)) for bp in ['A', 'G', 'C', 'T', 'N']}  for lst in product(['A', 'G', 'C', 'T', 'N'], repeat = 4)}
# cdef dict convdict1 = {ord('A'): 0, ord('C'): 2, ord('G'): 1, ord('T'): 3, ord('N'): 4}
fmerdict = {''.join(tmer):i for i, tmer in enumerate(product(['A', 'C', 'G', 'T', 'N'], repeat = 4))}
cdef int NNNN = fmerdict['NNNN']
cdef int REFBPN = 4
cdef int *fmer_ar = <int *> PyMem_Malloc(3125 * cython.sizeof(int))
if fmer_ar is NULL:
    raise MemoryError()
for tmer in product(['A', 'C', 'G', 'T', 'N'], repeat = 4):
    for b_i, b in enumerate(['A', 'G', 'C', 'T', 'N']):
        fmer_ar[fmerdict[''.join(tmer)] * 5 + b_i] = fmerdict[''.join(tmer[1:] + (b, ))]

cdef int *conv_arr = <int *> PyMem_Malloc(30 * cython.sizeof(int))
if conv_arr is NULL:
    raise MemoryError()
cdef int x = 0
while x < 30:
    conv_arr[x] = 100
    x += 1
cdef int A_C = ord('A')
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
conv_arr[0] = 0 #'A'
conv_arr[6] = 1 #'G'
conv_arr[2] = 2 #'C'
conv_arr[19] = 3 #'T'
conv_arr[13] = 4 #'N'

@cython.boundscheck(False)
cdef int conv_ar(int c):
    return conv_arr[c - A_C]
#     c = c - A_C
#     if c < 0 or c > 19:
#         raise AssertionError
#     cdef int v = conv_arr[c]
#     if v == 100:
#         raise AssertionError
#     return v

cpdef int validate_fmerdict(dict fmerdict):
    for k in fmerdict.keys():
        assert fmerdict[k] == fmer_ar[fmer_ar[fmer_ar[fmer_ar[NNNN * 5 + conv_ar(ord(k[0]))] * 5 + conv_ar(ord(k[1]))] * 5 + conv_ar(ord(k[2]))] * 5 + conv_ar(ord(k[3]))]
        for b in 'AGCTN':
            assert fmerdict[k] == fmer_ar[fmerdict[b + k[:-1]] * 5 + conv_ar(ord(k[-1]))]

cdef FTYPE_t THRESH = 5000

cpdef map_qual_score(char *seq, quals, char *gigar, string dest_id, int pos, int strand, int lseq, int mismatches, CSeqDict.CSeqDict dc,
                    np.ndarray[FTYPE_t, ndim=5] mdl0, np.ndarray[FTYPE_t, ndim=5] jnd0,
                    np.ndarray[FTYPE_t, ndim=4] mdl0s5, np.ndarray[FTYPE_t, ndim=4] jnd0s5,
                    np.ndarray[FTYPE_t, ndim=4] mdl0s2, np.ndarray[FTYPE_t, ndim=4] jnd0s2, 
                    np.ndarray[FTYPE_t, ndim=3] mdl0s25, np.ndarray[FTYPE_t, ndim=3] jnd0s25,
                    np.ndarray[FTYPE_t, ndim=3] mdl0s12, np.ndarray[FTYPE_t, ndim=3] jnd0s12,
                    np.ndarray[FTYPE_t, ndim=2] mdl0s125, np.ndarray[FTYPE_t, ndim=2] jnd0s125,
                    np.ndarray[FTYPE_t, ndim = 1] refins, double refins_lgsum, np.ndarray[FTYPE_t, ndim = 1] rdins, double rdins_lgsum):
    cdef char *curref
    cdef int ref_of, end, ofst, lga
    if strand == 1:
        curref = dc.get_seq(dest_id, 1)
    else:
        curref = dc.get_rev(dest_id)
    cdef long len_ref = dc.get_len(dest_id)
    
    if '>' in gigar:
        if strand == 1:
            ref_of = pos - 1
            end = pos - 1 + lseq
        else:
            ref_of = len_ref - pos + 1 - lseq 
            end = len_ref - pos + 1
         
        ga = gem_qa.get_opt_gigar_array_ptr(len_ref, curref, ref_of%len_ref, end%len_ref, strand, gigar, lseq, seq, mismatches, &lga, &ofst) 
    else:
        ofst = 0
        ga = gem_qa.create_gigar_array_ptr(gigar, &lga)
    if strand == 1:
        ref_of = pos - 1 + ofst
        end = pos - 1 + lseq + ofst + gem_qa.sum_skips2(ga, lga)
    else:
        ref_of = len_ref - pos  + 1 - lseq - ofst - gem_qa.sum_skips2(ga, lga)
        end = len_ref - pos + 1 - ofst 
    ref_of = ref_of % len_ref
    end = end % len_ref
    
    cdef int debug = 0
    cdef int prevmer = NNNN
    cdef double res = 0
    cdef int r_i = 0
    cdef int d_i = 0
    cdef int g_i = 0
    cdef int stupid = 1
    cdef int gb, inslen, convbp, i, cq, convrefbp
    cdef FTYPE_t curtot_sum, v
    
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
            prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
            r_i += 1
            while g_i+1 < lga and ga[g_i+1] == -1:
                g_i += 1
                inslen += 1
                prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
                r_i += 1
            res += rdins[inslen] - rdins_lgsum
        else:
            convbp = conv_arr[seq[r_i] - A_C]
            convrefbp = conv_arr[curref[ref_of+d_i] - A_C]
            if convrefbp == REFBPN:
                convrefbp = convbp
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
            
            res += log10(max_1(v)) - log10(curtot_sum)
            d_i += 1
            if ref_of + d_i >= len_ref:
                ref_of = (ref_of + d_i)%len_ref - d_i
            r_i += 1
            prevmer = fmer_ar[prevmer * 5 + convbp]
        g_i += 1
    PyMem_Free(ga)
    return 10**res

@cython.boundscheck(False)
cpdef double inner_func(unsigned char *seq, unsigned char *curref, np.ndarray[int, ndim = 1] ga, np.ndarray[FTYPE_t, ndim=4] tot, \
                        np.ndarray[FTYPE_t, ndim = 3] tots1, np.ndarray[FTYPE_t, ndim = 3] tots3, \
                        np.ndarray[FTYPE_t, ndim = 2] tots13, np.ndarray[FTYPE_t, ndim = 1] tots013, \
                        np.ndarray[FTYPE_t, ndim = 2] tots01, np.ndarray[FTYPE_t, ndim = 1] refins, \
                        double refins_lgsum, np.ndarray[FTYPE_t, ndim = 1] rdins, double rdins_lgsum):#, dict fmerdict):
    cdef int prevmer = NNNN
    cdef double res = 0
    cdef int r_i = 0
    cdef int d_i = 0
    cdef int g_i = 0
    cdef int lga = ga.shape[0]
    cdef int gb, inslen, convbp, i#, fmer
    cdef FTYPE_t curtot_sum, v
    
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
            prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
            r_i += 1
            while g_i+1 < lga and ga[g_i+1] == -1:
                g_i += 1
                inslen += 1
                prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
                r_i += 1
            res += rdins[inslen] - rdins_lgsum
        else:
            convbp = conv_arr[seq[r_i] - A_C]
            curtot_sum = tots3[r_i, prevmer, convbp]
            if curtot_sum < THRESH:
                curtot_sum = tots13[r_i, convbp]
                if curtot_sum < THRESH:
                    curtot_sum = tots013[convbp]
                    if curtot_sum < THRESH:
                        return 100
                    else:
                        v = tots01[convbp, conv_arr[curref[d_i] - A_C]]
                else:
                    v = tots1[r_i, convbp, conv_arr[curref[d_i] - A_C]]
            else:
                v = tot[r_i, prevmer, convbp, conv_arr[curref[d_i] - A_C]]
            res += log10(max_1(v)) - log10(curtot_sum)
            d_i += 1
            r_i += 1
            prevmer = fmer_ar[prevmer * 5 + convbp]
        g_i += 1
    return res

@cython.boundscheck(False)
cdef double inner_func2(char *seq, char *curref, int ref_of, int len_ref, int *ga, np.ndarray[FTYPE_t, ndim=4] tot, \
                        np.ndarray[FTYPE_t, ndim = 3] tots1, np.ndarray[FTYPE_t, ndim = 3] tots3, \
                        np.ndarray[FTYPE_t, ndim = 2] tots13, np.ndarray[FTYPE_t, ndim = 1] tots013, \
                        np.ndarray[FTYPE_t, ndim = 2] tots01, np.ndarray[FTYPE_t, ndim = 1] refins, \
                        double refins_lgsum, np.ndarray[FTYPE_t, ndim = 1] rdins, double rdins_lgsum, 
                        int lga, int debug) except 10000:#, dict fmerdict):
    cdef int prevmer = NNNN
    cdef double res = 0
    cdef int r_i = 0
    cdef int d_i = 0
    cdef int g_i = 0
    cdef int gb, inslen, convbp, i#, fmer
    cdef FTYPE_t curtot_sum, v
#     if debug == 1:
#         log_.warning(str(ref_of)); sleep(0.1)
    while g_i < lga:
        gb = ga[g_i]
#         if debug == 1:
#             log_.warning(str(g_i) + ': ' + str(gb))
        if gb == 1:
            inslen = 1
            while g_i+1 < lga and ga[g_i+1] == 1:
                g_i+=1
                inslen += 1
            res += refins[inslen] - refins_lgsum
            d_i += inslen
        elif gb == -1:
            inslen = 1
            prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
            r_i += 1
            while g_i+1 < lga and ga[g_i+1] == -1:
                g_i += 1
                inslen += 1
                prevmer = fmer_ar[prevmer * 5 + conv_arr[seq[r_i] - A_C]]
                r_i += 1
            res += rdins[inslen] - rdins_lgsum
        else:
#             if debug == 1:
#                 log_.warning(chr(seq[r_i])); sleep(0.1)
            convbp = conv_arr[seq[r_i] - A_C]
#             if debug == 1:
#                 log_.warning(convbp)
            curtot_sum = tots3[r_i, prevmer, convbp]
#             if debug == 1:
#                 log_.warning(curtot_sum)
            if curtot_sum < THRESH:
                curtot_sum = tots13[r_i, convbp]
#                 log_.warning(curtot_sum)
                if curtot_sum < THRESH:
                    curtot_sum = tots013[convbp]
#                     log_.warning(curtot_sum)
                    if curtot_sum < THRESH:
                        return 100
                    else:
                        v = tots01[convbp, conv_arr[curref[ref_of+d_i] - A_C]]
#                         log_.warning(v)
                else:
                    v = tots1[r_i, convbp, conv_arr[curref[ref_of+d_i] - A_C]]
#                     log_.warning(v)
            else:
#                 if debug == 1:
#                     log_.warning('before v %s %s %s %s %s' % (r_i, prevmer, convbp, ref_of+d_i, curref[ref_of+d_i]))#+ str(tot[r_i, prevmer, convbp, conv_arr[curref[ref_of+d_i] - A_C]]))
                v = tot[r_i, prevmer, convbp, conv_arr[curref[ref_of+d_i] - A_C]]
#                 if debug == 1:
#                     log_.warning('after v')
#                 log_.warning(v)
#             if debug == 1:
#                 log_.warning(str(v) + ' ' + str(curtot_sum)); sleep(0.1)
            res += log10(max_1(v)) - log10(curtot_sum)
            d_i += 1
            if ref_of + d_i >= len_ref:
#                 log_.warning('change %s %s' % (ref_of, ref_of + d_i))
                ref_of = (ref_of + d_i)%len_ref - d_i
#                 log_.warning('change %s %s' % (ref_of, ref_of + d_i))
            r_i += 1
            prevmer = fmer_ar[prevmer * 5 + convbp]
        g_i += 1
    return res

cdef char US = '_'
cdef char PLUS = '+'
cdef char * subst(char * st) except NULL:
    cdef char * pos = strrchr(st, US)
    cdef int ln = pos - st
    cdef char * res = <char *>PyMem_Malloc(ln + 1)
#     if res == NULL:
#         log_.warning('Memory in subst')
    strncpy(res, st, ln);
    res[ln] = '\0';
    return res
@cython.boundscheck(False)
cpdef dict do_everything(sriter, breakafter, np.ndarray[FTYPE_t, ndim=4] tot, \
                        np.ndarray[FTYPE_t, ndim = 3] tots1, np.ndarray[FTYPE_t, ndim = 3] tots3, \
                        np.ndarray[FTYPE_t, ndim = 2] tots13, np.ndarray[FTYPE_t, ndim = 1] tots013, \
                        np.ndarray[FTYPE_t, ndim = 2] tots01, np.ndarray[FTYPE_t, ndim = 1] refins, \
                        double refins_lgsum, np.ndarray[FTYPE_t, ndim = 1] rdins, double rdins_lgsum, 
                        CSeqDict.CSeqDict dc, dict dest_dict):
#     log_.warning('started'); sleep(0.1)
    cdef int c = 0
    cdef dict ret = {}
    cdef dict curret
    cdef char * cureid = 'not null'
    cdef char * eid 
    cdef int debug = 1 #$%$#
    for sr in sriter:
#         log_.warning('*****NEW SR' + sr.rid); sleep(0.1)
        eid = subst(sr.rid) #NEED TO FREE
        if strcmp(eid,cureid) != 0:
            curret = {}
#             PyMem_Free(cureid) #this will leak and I don't care
            cureid = subst(sr.rid)
            ret[cureid] = curret
#         log_.warning('sr' + sr.rid + ': do single sr'); sleep(0.1)
#         debug = int(sr.rid == 'srgca|GCA_000188735.1|asm|ASM18873v1|CP002489.1_0002')
        do_single_sr(sr.seq[0], len(sr.seq[0]), eid, curret, sr.se_maps, dc, dest_dict, tot,
                     tots1, tots3, tots13, tots013, tots01, refins, refins_lgsum, rdins, rdins_lgsum, debug)
        debug = 0
        c += 1
        if c % 10000 == 0:
            log_.info(c)
        if c > 10000 and breakafter:
#             log_.warning('C is'+str(c))
            return ret
        PyMem_Free(eid)
#     log_.warning('C is' + str(c))
    return ret
         
@cython.boundscheck(False)
cdef int do_single_sr(char *seq, int lseq, char *eid, dict curret, list se_maps, CSeqDict.CSeqDict dc,
                  dict dest_dict,#can fix? 
                  np.ndarray[FTYPE_t, ndim=4] tot, \
                        np.ndarray[FTYPE_t, ndim = 3] tots1, np.ndarray[FTYPE_t, ndim = 3] tots3, \
                        np.ndarray[FTYPE_t, ndim = 2] tots13, np.ndarray[FTYPE_t, ndim = 1] tots013, \
                        np.ndarray[FTYPE_t, ndim = 2] tots01, np.ndarray[FTYPE_t, ndim = 1] refins, \
                        double refins_lgsum, np.ndarray[FTYPE_t, ndim = 1] rdins, double rdins_lgsum, int debug) except -1:
    cdef int broken = 0, sawself = 0
    cdef dict curread = {}
    cdef char * dest_id
    cdef char * rseq
    cdef char * gigar
    cdef int strand, rlen, start, end, pos, ofst, len_ga
    cdef int *ga
    cdef double cur_p, max_other_p = -2000, curv, refv
    cdef char * dd_self = dest_dict[eid]
    cdef char * dd_other
    for m in se_maps:
#         debug = int(m.dest_id == 'gca|GCA_000211395.2|asm|ASM21139v1|AFET01000004.1')
#         if debug == 1:
#             log_.warning('%s %s %s %s' % (m.dest_id, m.gigar, m.pos, m.strand)); sleep(0.1)
        dest_id = m.dest_id
        gigar = m.gigar
        pos = m.pos
        if strcmp(dest_id, eid) == 0:
            sawself = 1
        if m.strand == u'+':
            strand = 1
            rseq = dc.get_seq(m.dest_id, 1)
        else:
            strand = -1
            rseq = dc.get_rev(m.dest_id)
        rlen = dc.get_len(m.dest_id)
        if '>' in m.gigar:
            if strand == 1:
                start = pos - 1
                end = pos - 1 + lseq
            else:
                start = rlen - pos + 1 - lseq 
                end = rlen - pos + 1
             
            ga = gem_qa.get_opt_gigar_array_ptr(rlen, rseq, start%rlen, end%rlen, strand, m.gigar, lseq, seq, m.mismatches, &len_ga, &ofst) 
        else:
            ofst = 0
            ga = gem_qa.create_gigar_array_ptr(gigar, &len_ga)
        if strand == 1:
            start = pos - 1 + ofst
            end = pos - 1 + lseq + ofst + gem_qa.sum_skips2(ga, len_ga)
        else:
            start = rlen - pos  + 1 - lseq - ofst - gem_qa.sum_skips2(ga, len_ga)
            end = rlen - pos + 1 - ofst 
        start = start % rlen
        end = end % rlen
        
#         if debug == 1:
# #             log_.warning('inner_func2'); sleep(0.1)
#             if m.dest_id == 'gca|GCA_000211395.2|asm|ASM21139v1|AFET01000004.1':
#                 log_.warning('now!'); sleep(0.1)
        cur_p = inner_func2(seq, rseq, start, rlen, 
                                                 ga, tot, tots1, tots3, \
                                                 tots13, tots013, tots01, refins, refins_lgsum, \
                                                 rdins, rdins_lgsum, len_ga, 0)
#         if debug == 1:
#             log_.warning('inner_func2:done'); sleep(0.1)
        if cur_p == 100:
            broken = 1
            break
        dd_other = dest_dict[dest_id]
        if dd_other != dd_self and cur_p > max_other_p:
            max_other_p = cur_p
        if dd_other not in curread or cur_p > curread[dd_other]:
            curread[dd_other] = cur_p
        PyMem_Free(ga)
    
    if broken == 1 or sawself == 0 or max_other_p > curread[dd_self]:
        return 0
    refv = curread[dd_self]
    for k, v in curread.iter():
        curv = v
        curv = pow(10, curv - refv)
        curret[k] = curret.get(k, 0) + curv
    return 0
        
