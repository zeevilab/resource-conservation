from os.path import join, exists, basename
from glob import glob
from Data.config import SNP, Calling
from Data.SNP.SNPAnalysis import getvariants, analyze_genes
from lib import Utils
from Data.SNP.CollatedGeneGroups import filter_db, do_one_group_pnps, do_collate
from lib.Utils import split3way

def main():
    # Filter VCF files for high quality SNPs and save them in DataFrame
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    # estimated total CPU time for this part ~1,000 hours (Intel(R) Xeon(R) CPU E5-2690 v3)
    for fname in glob(join(SNP.OM_RGC.InputDir, '*.vcf')):
        if exists(join(SNP.OM_RGC.GeneDFDir, basename(fname).replace('.vcf','.df'))):
            continue
        getvariants(fname, SNP.OM_RGC.GeneDFDir,only_snps=SNP.OM_RGC.OnlySNPs, 
                    qual_thresh=SNP.OM_RGC.QualThresh, min_samples=SNP.OM_RGC.MinSamples, 
                    min_variants=SNP.OM_RGC.MinVariants)
    
    # Calculate selection metrics per gene (e.g. pN/pS)
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    # estimated total CPU time for this part >5,000 hours
    for fname in glob(join(SNP.OM_RGC.GeneDFDir, '*.df')):
        outdir = SNP.OM_RGC.OutDir
        if all([exists(join(outdir, basename(fname).replace('.df',ext)))\
                for ext in ['.pnps.df', '.ffdeg_pi_wit.df', '.pnpn.df']]):
                    continue
        analyze_genes(fname, Calling.OM_RGC.DbFasta, outdir, SNP.OM_RGC.CacheDir,
                      min_pos_reads=SNP.OM_RGC.MinPosReads, 
                      min_perc_poss=SNP.OM_RGC.MinPosReads, 
                      min_total_var_support=SNP.OM_RGC.MinTotalVarSupport, 
                      min_maf=SNP.OM_RGC.MinMaf,
                      min_samples=SNP.OM_RGC.MinSamples, 
                      min_variants=SNP.OM_RGC.MinVariants)
    
    # Collate selection metrics per KEGG KO / eggNOG OG. Requires running eggnogMapper on OM-RGC
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    # Note: this calculates only pN/pS per KEGG KO / eggNOG OG. To calculate fourfold degenerate
    # pi within (validation) or pN(conservative AA substitutions) vs pN(radical AA substitutions)
    # (also validation), refer to the relevant methods within SNP/CollatedGeneGroups.py
    for db in SNP.OM_RGC.GeneGroupCollateDBs:
        dbdct = Utils.Load(db)
        dbdct = filter_db(dbdct, SNP.OM_RGC, SNP.OM_RGC.MinGenes)
        for nm, genes in dbdct.items():
            if not exists(join(SNP.OM_RGC.OutDirCollate, 'pnps',
                               split3way(db)[1] + '_' + nm + '.pnps.df')):
                do_one_group_pnps(nm, split3way(db)[1], genes,SNP.OM_RGC, SNP.OM_RGC.MinGenes)
                
    # This groups all of the files and saves them in the output folder defined as General.Basepath
    do_collate(join(SNP.OM_RGC.OutDirCollate, 'pnps', 'KEGG'), 4, 60, 5, 50, 5)
    do_collate(join(SNP.OM_RGC.OutDirCollate, 'pnps', 'eggNOG'), 4, 60, 5, 50, 5)

if __name__ == '__main__':
    main()