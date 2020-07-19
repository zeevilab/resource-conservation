from collections import namedtuple
from operator import attrgetter
import logging
import cy_ext_seq.gem_qa as gqa
import ujson
from gzip import open as opengz
from tools.sam2pmp import SourceReadSAM
log_ = logging.getLogger(__name__)

PE_DELIM = '::'
ALIGN_DELIM = ','
MATCHBLOCK_DELIM = ':'
PLUS = '+'
MINUS = '-'

DestMap_se = namedtuple('DestMap_se', ['dest_id', 'mismatches', 'strand', 'pos', 'gigar', 'end'])
DestMap_pe = namedtuple('DestMap_pe', ['dest_id', 'mismatches1', 'mismatches2', 'strand1', 'strand2', 'pos1', 'pos2', 'gigar1', 'gigar2'])

class SourceRead(object):
    def __init__(self, rid, seq, qual):
        self.rid = rid
        self.seq = [str(s) for s in seq]  # #TODO: remove
        self.quals = qual
        self.se_maps = []
        self.pe_maps = []
    
    def to_ser(self):
        return (self.rid, self.seq, self.quals, [tuple(x) for x in self.se_maps], [tuple(y) for y in self.pe_maps])
    
    @staticmethod
    def from_ser(ser):
        res = SourceRead(ser[0], ser[1], ser[2])
        res.se_maps = [DestMap_se(*x) for x in ser[3]]
        res.pe_maps = [DestMap_pe(*y) for y in ser[4]]
        return res
    
    def __len__(self):
        return len(self.se_maps) + len(self.pe_maps)
    
    def sort(self):
            self.se_maps.sort(key=attrgetter('mismatches'))
            self.pe_maps.sort(key=lambda m: m.mismatches1 + m.mismatches2)

    #
    def try_fix_pos(self, gd, seq, gid, strand, pos, gigar):
        for off in (1, -1, 2, -2, 3, -3):
            w = _is_wrong_strand(gd, seq, gid, strand, pos + off, gigar)
            if w != 2:
                return w, off
        return 2, 0
    #  
    def _add_from_cells(self, cells, end, seqdict):
        assert len(cells) == 5 and cells[0].startswith(self.rid)
        if cells[4] == MINUS:
            return
        alignments = cells[4].split(ALIGN_DELIM)
        for align in alignments:
            if '\x00' in align:
                log_.warning('Encountered null character in alignment ' + str(alignments.index(align)) + ' for read ' + cells[0] + '. Skipping!: ' + repr(align))
                continue
            if PE_DELIM in align:
                tmpalgn = align.split(PE_DELIM)
                tmpalgn[0] = tmpalgn[0].split(MATCHBLOCK_DELIM, 4)
                tmpalgn[1] = tmpalgn[1].split(MATCHBLOCK_DELIM, 4)
                ga1 = gqa.create_gigar_array(tmpalgn[0][3])
                ga2 = gqa.create_gigar_array(tmpalgn[1][3])
                if tmpalgn[0][1] == tmpalgn[1][1] or tmpalgn[0][0] != tmpalgn[1][0]:
                    seqs = cells[1].split()
                    w1 = _is_wrong_strand(seqdict, seqs[0], tmpalgn[0][0].split()[0], tmpalgn[0][1], int(tmpalgn[0][2]), ga1)
                    w2 = _is_wrong_strand(seqdict, seqs[1], tmpalgn[1][0].split()[0], tmpalgn[1][1], int(tmpalgn[1][2]), ga2)
                    if w1 == 1 and w2 == 0:
                        tmpalgn[0][1] = PLUS if tmpalgn[0][1] == MINUS else MINUS
                    elif w2 == 1 and w1 == 0:
                        tmpalgn[1][1] = PLUS if tmpalgn[1][1] == MINUS else MINUS
                    elif w1 != 0 or w2 != 0:
                        if w1 == 2:
                            w1, off = self.try_fix_pos(seqdict, seqs[0], tmpalgn[0][0].split()[0], tmpalgn[0][1], int(tmpalgn[0][2]), ga1)
                            if w1 != 2:
                                tmpalgn[0][2] = str(int(tmpalgn[0][2]) + off)
                                log_.info('Fixed a read by offest ' + str(off) + '\n\t' + seqs[0] + '\n\t' + seqs[1] + '\n\t' + align)
                        if w2 == 2:
                            w2, off = self.try_fix_pos(seqdict, seqs[1], tmpalgn[1][0].split()[0], tmpalgn[1][1], int(tmpalgn[1][2]), ga2)
                            if w2 != 2:
                                tmpalgn[1][2] = str(int(tmpalgn[1][2]) + off)
                                log_.info('Fixed a read by offest ' + str(off) + '\n\t' + seqs[0] + '\n\t' + seqs[1] + '\n\t' + align)
                        if w1 == 1 and w2 == 0:
                            tmpalgn[0][1] = PLUS if tmpalgn[0][1] == MINUS else MINUS
                        elif w2 == 1 and w1 == 0:
                            tmpalgn[1][1] = PLUS if tmpalgn[1][1] == MINUS else MINUS
                        elif w1 != 0 or w2 != 0:  
                            log_.info('This pe read is broken' + '\n\t' + seqs[0] + '\n\t' + seqs[1] + '.\n Skipping. This is probably not too bad unless you see too many of these.')
                if tmpalgn[0][0] != tmpalgn[1][0]:
                    log_.warning('Encountered paired mapping to different bacteria!')
                    log_.warning(cells)
                    self.se_maps.append(_se_from_algn(tmpalgn[0], 0, _get_mms_ga(ga1), seqdict))
                    self.se_maps.append(_se_from_algn(tmpalgn[1], 1, _get_mms_ga(ga2), seqdict))
                else:
                    self.pe_maps.append(_pe_from_algn(tmpalgn, seqdict, ga1, ga2))
            else:
                if end == 2:
                    raise AssertionError
                algn = align.split(MATCHBLOCK_DELIM, 4)
                ga = gqa.create_gigar_array(algn[3])
                self.se_maps.append(_se_from_algn(algn, end, _get_mms_ga(ga), seqdict))

def _is_wrong_strand(gd, seq, gid, strand, pos, gigar):
    matched_cor = _validatematch(gd, seq, gid, strand, pos, gigar)
    if matched_cor:
        return 0
    matched_wrong = _validatematch(gd, seq, gid, PLUS if strand == MINUS else MINUS, pos, gigar)
    if matched_wrong:
        return 1
    else:
        return 2
    
def _validatematch(gd, seq, gid, strand, pos, gigar_arr):
    lngt = len(seq) + gqa.sum_skips(gigar_arr)
    genesq, start, end, strand = gd.get_sequence(gid, strand, pos, lngt)
    destseq = _circ_extract_subseq(genesq, start, end)
    return gqa.validate_match3(str(seq), str(destseq), gigar_arr)

def _circ_extract_subseq(seq, start, end):
        if start < end:
            return seq[start:end]
        else:
            return seq[start:] + seq[:end]
        
def _se_from_algn(algn, end, gigar_arr, seqdict):
    destid = algn[0].split()[0]
    return DestMap_se(destid, gigar_arr, algn[1], (int(algn[2]) - 1) % seqdict.get_len(destid) + 1, algn[3], end)

def _pe_from_algn(algn, seqdict, ga1, ga2):
    destid = algn[0][0].split()[0]
    destlen = seqdict.get_len(destid)
    return DestMap_pe(destid, _get_mms_ga(ga1), \
                      _get_mms_ga(ga2), algn[0][1], algn[1][1], \
                      (int(algn[0][2]) - 1) % destlen + 1, \
                      (int(algn[1][2]) - 1) % destlen + 1, \
                      algn[0][3], algn[1][3])

def _get_mms_ga(g_a):
    return g_a.nonzero()[0].shape[0]

def _load_from_file(fpath, sam_based):
    with opengz(fpath) as fin:
        for r in _load_iterable_fromdesc(fin, sam_based):
            yield r

def _load_iterable_fromdesc(desc, sam_based):
    sr = SourceReadSAM if sam_based else SourceRead
    try:
        while True:
            yield sr.from_ser(ujson.loads(desc.readline()))
    except ValueError as ve:
        if ve.args[0] == 'No JSON object could be decoded' or \
            ve.args[0] == 'Expected object or value':
            raise StopIteration
        else:
            raise ve
    except EOFError:
        raise StopIteration
