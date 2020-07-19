import sys
import inspect
from os.path import join, basename
import logging
import numpy as np
from datetime import datetime
from collections import defaultdict
from pandas import to_pickle
import Data.Mapping.cy_ext_seq.gem_qa as gqa
import ujson
import gzip
from Data.Mapping.ReadContainer import ReadContainer, _get_flank_size
from lib.Utils import open_gz_indif
from Data.Mapping.pmap import _load_from_file
from Data.config import Mapping
log_ = logging.getLogger('ICRA')

OCEAN_GENES = 'oceangenes'
 
def single_file(fq1, fq2, outfol, max_mismatch, consider_lengths, epsilon, \
                max_iterations, min_bins = 4, max_bins = 100, min_reads = 5, 
                dense_region_coverage = 60, length_minimum = 300, \
                length_maximum = 1e6, usage = OCEAN_GENES, use_theta = False, pmpf = None, 
                average_read_length = None, force_save_delta = False, sam_based=False):
    
    length_db_f = Mapping.OM_RGC.LengthsFile if usage == OCEAN_GENES else None
    
    read_container, pi, theta1, average_read_length, lengthdb = \
        _initialize(fq1, fq2, pmpf, max_mismatch, consider_lengths, \
                    length_db_f, length_minimum, length_maximum, min_reads, min_bins, max_bins, \
                    dense_region_coverage, use_theta, average_read_length, sam_based)
    if len(pi) == 0:
        delta = {}
    else:
        delta, pi = _runIterative(pi, theta1, read_container, min_bins, max_bins, \
                              min_reads, average_read_length, max_mismatch, \
                              lengthdb, dense_region_coverage, consider_lengths, \
                              length_minimum, length_maximum, epsilon, max_iterations, use_theta)
    with gzip.open(join(outfol, basename(pmpf).replace('.pmp', '.jspi')), 'wt') as of:
        ujson.dump(pi, of, double_precision=100)
    with gzip.open(join(outfol, basename(pmpf).replace('.pmp', '.jsdel')), 'wt') as of:
        ujson.dump(delta, of, double_precision=100)
    return pmpf.replace('.pmp', '.jspi'), pmpf.replace('.pmp', '.jsdel')

def _initialize(fq1, fq2, pmpf, max_mismatch, consider_lengths, \
               length_db_f, length_minimum, length_maximum, min_reads, min_bins,  
               max_bins, dense_region_coverage, use_theta, average_read_length, 
               sam_based): 
    if average_read_length is None:
        average_read_length = _getaveragereadlength(fq1)
    lengthdb = _get_len_dct(None, None, length_db_f)
    read_container = ReadContainer()
    pe = fq2 is not None 
    allowed_dests = {k for k,v in lengthdb.items() if length_minimum <= v <= length_maximum}
    for sr in _load_from_file(pmpf, sam_based): 
        read_container.add_pmp_sr(sr, pe, allowed_dests, max_mismatch, None, None, None)
    del allowed_dests
    pi = read_container.get_pi()
    
    read_container.remove_dests([k for k,val in pi.items() if val<=min_reads])
    pi = {k:val for k, val in pi.items() if  val > min_reads}
    
    read_container.init_cov_dict(pi, min_bins, max_bins, min_reads, average_read_length, lengthdb)
    #note: this will create a fully updated coverage dict in the process
    pi = read_container.get_dense_region_coverage_pi(dense_region_coverage) 
    
    read_container.remove_dests([k for k,val in pi.items() if val<=min_reads])
    pi =  {k:v for k,v in pi.items() if v>min_reads}
    
    if consider_lengths:
        flank_sizes = _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length)
        pi = {k:(v/(lengthdb[k] - (2 * flank_sizes[k]))) for k, v in pi.items()}
    pisum = sum(pi.values())
    pi = {k:v/float(pisum) for k,v in pi.items()}
    
    theta1 = read_container.get_theta() if use_theta else None
    
    return read_container, pi, theta1, average_read_length, lengthdb
    
def _runIterative(pi, theta1, read_container, min_bins, max_bins, min_reads, average_read_length, \
                  max_mismatch, lengthdb, dense_region_coverage, consider_lengths, length_minimum, 
                  length_maximum, epsilon, max_iterations, use_theta):
    starttime = datetime.now()
    pi_dict = {}
    i = 1
    while True:
        pi_dict[i] = pi
        delta = read_container.do_delta(pi, theta1, 1. / ((len(pi)**2) * (10**max_mismatch)), use_theta)
        read_container.recalc_cov_bins()
        prevPi = pi
        #reminder: implicitly creates an updated covdic
        pi = read_container.get_dense_region_coverage_pi(dense_region_coverage, delta)
        read_container.remove_dests([k for k,val in pi.items() if val<min_reads])
        pi =  {k:v for k,v in pi.items() if v>=min_reads}
        
        if len(pi) == 0:
            log_.info("No adequately covered strains found")
            return {}, {}
        
        theta1 = read_container.get_theta() if use_theta else None
        
        if consider_lengths:
            flank_sizes = _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length)
            pi = {k:(v/(lengthdb[k] - (2 * flank_sizes[k]))) for k, v in pi.items()}
        pisum = sum(pi.values())
        pi = {k:v/float(pisum) for k,v in pi.items()}
        
        prevPi = {k:v for k, v in prevPi.items() if k in pi}
        dPi = _LogDistDict(pi, prevPi)
         
        i += 1
        log_.info("Iteration {} - Time: {}, dPi = {:.2e}, nPi = {}".format(i, datetime.now() - starttime, dPi, 
                                                                            len(pi)))
         
        if  dPi < epsilon or i > max_iterations:
            break
    log_.info("Final result - Time: {}".format(datetime.now() - starttime))
    return read_container.get_full_delta(pi, theta1, 1. / ((len(pi)**2) * (10**max_mismatch)), use_theta), pi

def _write_if_vb(vb, pth, obj):
    if vb:
        to_pickle(obj, pth)

def _logdist(u, v):
    x = np.abs(np.log10(np.asarray(u)) - np.log10(np.asarray(v)))
    return np.percentile(x, 90)

def _LogDistDict(u, v):
    x = []
    y = []
    if not (set(u.keys()).intersection(set(v.keys())) == set(u.keys()) == set(v.keys())):
        raise KeyError("Keys not overlapping between the two dictionaries")
    for key in u:
        if u[key]>1e-5 or v[key]>1e-5:
            x.append(u[key])
            y.append(v[key])
    return _logdist(x,y)

def _getaveragereadlength(fq):
    with open_gz_indif(fq) as inf:
        sum_bps = 0
        for i, read in enumerate(inf):
            if i % 4 == 1:
                sum_bps += len(read) - 1 #has a /n at the end
            if i >= 4000: break
    arl = int(np.round(sum_bps/1000))
    log_.info("Average read length for {} set to {}".format(basename(fq), arl))
    return arl

def _get_len_dct(sequence_dict = None, read_dest_dict = None, length_db_f = None):
    if length_db_f is not None:
        with open(length_db_f, 'rt') as fin:
            lengthdb =  ujson.load(fin)
    elif hasattr(sequence_dict, 'get_lendct'):
        lengthdb = sequence_dict.get_lendct()
    else:
        lengthdb = {k:sequence_dict.get_len(k) for k in sequence_dict._dct.keys()}
    if read_dest_dict is not None:
        newlengthdb = defaultdict(int)
        for l in lengthdb:
            newlengthdb[read_dest_dict[l][0]] += lengthdb[l]
        return dict(newlengthdb)
    else:
        return lengthdb

def _get_opt_gigar_array(seqdict, dest_id, strand, pos, mismatches, seq, gigar):
    if '>' in gigar:
        lseq = len(seq)
        genesq, start, end, strand = seqdict.get_sequence(dest_id, strand, pos, lseq)
        return gqa.get_opt_gigar_array(seqdict.get_len(dest_id),genesq, start, end,
                                                                 strand, gigar, lseq, seq, mismatches)
    else:
        return gqa.create_gigar_array(gigar), 0
  
def get_size(obj, seen=None):
    """Recursively finds size of objects in bytes"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if hasattr(obj, '__dict__'):
        for cls in obj.__class__.__mro__:
            if '__dict__' in cls.__dict__:
                d = cls.__dict__['__dict__']
                if inspect.isgetsetdescriptor(d) or inspect.ismemberdescriptor(d):
                    size += get_size(obj.__dict__, seen)
                break
    if isinstance(obj, dict):
        size += sum((get_size(v, seen) for v in obj.values()))
        size += sum((get_size(k, seen) for k in obj.keys()))
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum((get_size(i, seen) for i in obj))
    return size

def _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length):
    numbins = {k:max(min_bins, min(max_bins,int(v/100))) + 2 for k,v in pi.items()}
    return {k:_get_flank_size(lengthdb[k], numbins[k], average_read_length) for k,v in pi.items()}

