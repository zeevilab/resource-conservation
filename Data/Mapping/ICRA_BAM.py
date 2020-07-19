from gzip import open as gzopen
import ujson
import pysam
import logging
from _datetime import datetime
from itertools import groupby
from operator import attrgetter
log_ = logging.getLogger('ICRA_BAM')

def _write_to_sam(sam_fh, alngrp):
    for aln in alngrp:
        sam_fh.write(aln)
        
def _add_zw_tag(aln, val):
    tags = aln.get_tags(with_value_type=True)
    tags.append(('ZW',val,'f'))
    aln.set_tags(tags)
    return aln

def _add_deltas(alngrp, d_rid, delta_thresh):
    d_rid_dct = {}
    for d in d_rid:
        d_rid_dct[(d[0],d[1])] = d[3]
        if d[2] > 0: #PE
            d_rid_dct[(d[0],d[2])] = d[3]
    for aln in alngrp:
        try:
            prob = d_rid_dct[(aln.reference_name, aln.pos)]
        except KeyError:
            prob = 0
        prob = 1 if prob > delta_thresh else 0 if prob < (1-delta_thresh) else prob
        aln = _add_zw_tag(aln, prob)
    return alngrp
        
def add_ICRA_probs(jsdel_f, in_bam_f, out_bam_f, remove_unmapped=True, remove_not_in_delta=False,
                      delta_thresh=0.9):
    dt = datetime.now()
    with gzopen(jsdel_f, 'rt') as jsdel_fh:
        delta = dict(ujson.load(jsdel_fh, precise_float = True))
    log_.info('Loaded delta from file {}. Time: {}'.format(jsdel_f, datetime.now()-dt))
    in_bam = pysam.AlignmentFile(in_bam_f)  # @UndefinedVariable
    out_bam = pysam.AlignmentFile(out_bam_f, "wb", header=in_bam.header)  # @UndefinedVariable
    for rid, grp in groupby(in_bam, attrgetter('query_name')):
        alngrp = list(grp)
        if remove_unmapped:
            if (len(alngrp) == 1 and alngrp[0].is_unmapped) \
                    or (len(alngrp) == 2 and alngrp[0].is_unmapped and alngrp[1].is_unmapped):
                continue
        try:
            d_rid = delta[rid]
        except KeyError:
            if remove_not_in_delta:
                continue
            for aln in alngrp:
                aln = _add_zw_tag(aln, 0)
            _write_to_sam(out_bam, alngrp)
            continue
        alngrp = _add_deltas(alngrp, d_rid, delta_thresh)
        _write_to_sam(out_bam, alngrp)
    
