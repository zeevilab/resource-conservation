from glob import glob
from os.path import join, basename, exists, isdir
import os
import logging
import pysam
from pandas.io.pickle import read_pickle
from pandas.core.reshape.concat import concat
from lib.Utils import shell_command, rmrf, mkdirifnotexists
from Data.config import Calling, Biodata
from Data.Calling.FilterICRABAM import do
from Data.Calling.SortMpileupCall import mpilupcall, samcat, samsort, splitbam
log_ = logging.getLogger('CallingOcean')
del logging

THREADS = 4

def filt_sort_split(fnames, outfol, out_fname, reference_list, prob_thresh, threads):
    filts = []
    unite_f = join(outfol, out_fname + '.u.bam')
    sort_f = unite_f.replace('.u.bam','.s.bam')
    if not exists('{}.done'.format(join(outfol, out_fname))):
        for fname in fnames:
            filt_f = join(outfol, basename(fname.replace('.icra.bam','.filt.bam')))
            filts.append(filt_f)
            do(fname, filt_f, prob_thresh)
        samcat(filts, unite_f, delinfile=True)         
        samsort(unite_f, sort_f, threads, delinfile=True)
        shell_command('touch {}.done'.format(join(outfol, out_fname)))
    splitbam(sort_f, outfol, reference_list)

def mpileup_call(dirnm, threads, fasta_db, mpileup_params, del_indir=True):
    shell_command('ulimit -n 16384')
    outvcf = dirnm + '.vcf'
    donef = dirnm + '_calling.done'
    if exists(donef):
        return
    mpilupcall([join(dirnm, '*.bam')], fasta_db, outvcf, 
               Q=mpileup_params['min_base_qual'], 
               L=mpileup_params['max_depth_indel'], d=mpileup_params['max_depth'], 
               m=mpileup_params['min_iReads'], threads=threads)
    if del_indir:
        rmrf(dirnm)
    shell_command('touch {}'.format(donef))
        
def do_references(allbams, outfol, reference_list, prob_thresh, fasta_db, 
                  mpileup_params, threads):
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    for unite_fname, grp in allbams.groupby(level=0):
        fnames = grp.stack().values
        if all([exists(fname.replace('.icra.bam','.icrabamdone')) for fname in fnames]):
            filt_sort_split(fnames, outfol, unite_fname, reference_list, prob_thresh, threads-1)
        else:
            raise RuntimeError("Make sure mapping part is finished before running variant calling")
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline
    # Make sure you wait for jobs from previous loop to finish before you start this one
    for dirnm in sorted(glob(join(outfol, '*'))):
        if isdir(dirnm) and  dirnm.split('/')[-1] in reference_list: 
            mpileup_call(dirnm, threads-1, fasta_db, mpileup_params, True)

def main():
    os.chdir(mkdirifnotexists(join(Calling.OM_RGC.CallDir,'tmp')))
    bioG_m = read_pickle(Biodata.bioGEOTRACES.metadataDF)
    ALOHA_m = read_pickle(Biodata.ALOHA_BATS.metadataDF)
    TARA_m = read_pickle(Biodata.TARA.metadataDF)
    allbams = concat([TARA_m, ALOHA_m, bioG_m], sort=False)[['ICRABAM_1','ICRABAM_2']]
    dirnames = sorted(list(set([ref[:-4] for ref \
                        in pysam.AlignmentFile(allbams.iloc[-1]['ICRABAM_1']).header.references])))
    # This is set to process 80 genes (each with all samples) at a time.
    # Changing it to a higher setting will cause everything to run faster on a hpc system, but
    # take up more memory and space for intermediate files
    dirnamegrps = [dirnames[i:i+80] for i in range(0, len(dirnames), 80)]
    for reference_list in dirnamegrps:
        # TODO: IMPORTANT! Wrap the loops in the called method with your hpc job submission pipeline 
        # Also IMPORTANT! Make sure each loop runs synchronously with the next (wait for one to 
        # finish before you start the next)
        # Estimated total CPU time for this part >25,000 hours (Intel(R) Xeon(R) CPU E5-2690 v3)
        do_references(allbams, Calling.OM_RGC.CallDir, reference_list, 
                      Calling.OM_RGC.FilterThreshold, Calling.OM_RGC.DbFasta, 
                      Calling.OM_RGC.mpileupParams, THREADS)
            
if __name__ == '__main__':
    main()