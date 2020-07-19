from os.path import join, basename, exists
import os
from glob import glob
from lib.Utils import mkdirifnotexists
from Data.config import Mapping, RawFastq
from Data.Mapping.MapPmpICRA import do_single_mappmp, do_single_icrabam

def _create_length_db():
    if exists(Mapping.OM_RGC.LengthsFile):
        return
    lens = {}
    from Bio import SeqIO
    import ujson
    for rec in SeqIO.parse(Mapping.OM_RGC.IndexFasta, 'fasta'): 
        lens[rec.id] = len(rec.seq)
    with open(Mapping.OM_RGC.LengthsFile, 'wt') as fout:
        ujson.dump(lens, fout)
    del ujson
    del SeqIO
         
def map_sample(force_rerun=False):
    #Single run operation to create a lengths database
    _create_length_db()
    indexf=Mapping.OM_RGC.IndexFile
    os.chdir(mkdirifnotexists(join(Mapping.OM_RGC.MapDir, 'tmp')))
    all_files = glob(join(RawFastq.ALOHA_BATS.FastqDir, '*.fastq.gz')) + \
                glob(join(RawFastq.TARA.FastqDir, '*.fastq.gz')) + \
                glob(join(RawFastq.bioGEOTRACES.FastqDir, '*.fastq.gz'))
    flist = []
    for fq in all_files:
        prefix = basename(fq.replace('.fastq.gz',''))  
        prefix = join(Mapping.OM_RGC.MapDir, prefix)
        flist.append((fq, prefix))
        
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    # estimated total CPU time for this part >10,000 hours (Intel(R) Xeon(R) CPU E5-2690 v3)
    for fq, prefix in flist:
        if exists(prefix + '.icrabamdone') and not force_rerun:
            print('{} pipeline completed, no need to rerun'.format(prefix))
            continue
        args = (fq, None, prefix, indexf)
        mapkwargs = dict(threads=16, map_param_dict=Mapping.OM_RGC.MapParams)
        if not exists(prefix + '.pmpdone') or force_rerun:
            do_single_mappmp(*args, **mapkwargs)

    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    # estimated total CPU time for this part >25,000 hours
    for fq, prefix in flist:
        args = (fq, None, prefix, indexf)
        icrakwargs = dict(icra_usage=Mapping.OM_RGC.ICRAUsage, 
                          icra_param_dict=Mapping.OM_RGC.ICRAParams,
                          remove_unmapped=Mapping.OM_RGC.RemoveUnmapped, 
                          remove_not_delta=Mapping.OM_RGC.RemoveNotDelta, 
                          delta_thresh=Mapping.OM_RGC.DeltaThresh,
                          delete_pmp=Mapping.OM_RGC.DeletePMP, 
                          delete_old_bam=Mapping.OM_RGC.DeleteOldBAM)
        if exists(prefix + '.icrabamdone') and not force_rerun:
            continue
        # NOTE: this process may be memory intensive, especially for bam files larger than 25GB.
        # To avoid problems, process these files only on nodes with > 256GB memory
        do_single_icrabam(args, icrakwargs)

if __name__ == '__main__':
    map_sample()