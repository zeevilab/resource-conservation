from os.path import join, basename
import pysam
from itertools import groupby
from lib.Utils import shell_command, tryrm, mkdirifnotexists

def samcat(fnames, unite_f, delinfile=True):
    shell_command('samtools cat -o {} {}'.format(unite_f, ' '.join(fnames)))
    if delinfile:
        for fname in fnames:
            tryrm(fname)

def samsort(fname_in, fname_out, threads, delinfile = False):
    shell_command('samtools sort -@ {} -o {} {}'.format(threads-1, fname_out, fname_in))
    if delinfile:
        tryrm(fname_in)
    
def mpilupcall(fnames, db_fasta, outvcf, Q=15, L=1000, d=100000, m=2, threads=8):
    if len(fnames) == 0:
        return
    cmd = 'bcftools mpileup --threads {} -a FORMAT/AD -Q {} -L {} -d {} -m {} -f {} {}'\
            .format(threads, Q, L , d, m, db_fasta, ' '.join(fnames))
    cmd += ' | bcftools call --threads {} -Ov -mv -o {}'.format(threads, outvcf)
    shell_command(cmd)
    
def _dirfromrec(rec):
    return rec.reference_name[:-4]

def splitbam(bam_fname, outbasedir, reference_list):
    ## Counting on input to be sorted
    with pysam.AlignmentFile(bam_fname) as af_in:  # @UndefinedVariable
        for dirnm, grp in groupby(af_in, _dirfromrec):
            if reference_list is not None and dirnm not in reference_list:
                continue
            lgrp = list(grp)
            if len(lgrp) == 0: continue
            reference_names = [af_in.get_reference_name(i) for i in range(lgrp[-1].reference_id+1)]
            reference_lengths = [af_in.header.get_reference_length(nm) for nm in reference_names]
            outsplitdir = mkdirifnotexists(join(outbasedir, dirnm))
            with pysam.AlignmentFile(join(outsplitdir, basename(bam_fname.replace('.s.filt',''))),  # @UndefinedVariable
                                     'wb', reference_names=reference_names, 
                                     reference_lengths=reference_lengths) as af_out:
                for rec in lgrp:
                    af_out.write(rec)
