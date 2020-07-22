from os.path import join
from lib.Utils import mkdirifnotexists
from Data.Mapping.Bowtie2WrapperSlim import MapPreset
from Data.Mapping.ICRA import OCEAN_GENES

class General:
    #TODO: replace with filesystem base (the project requires ~35TB of storage for all intermediates
    Scratch = '~' 
    #TODO: replace with output path for finished dataframes
    Basepath = mkdirifnotexists(join(Scratch, 'Analyses/2019-Oceans/DFOut')) 
    #TODO: replace with output path for intermediates
    Tmppath = mkdirifnotexists(join(Scratch, 'Analyses/2019-Oceans/tmp'))
    #TODO: replace to where you place your data, i.e., samples and databases
    Data = join(Scratch, 'Data')
    #TODO: replace with where you place your samples (e.g. Tara)
    Samples = join(Data, 'Samples')
    #TODO: replace with where you place your sample metadata 
    Metadata = join(Data, 'Metadata')
    #TODO: replace with where you place your databases (e.g. OM-RGC) 
    Databases = join(Data, 'Databases')
    
    
class Biodata:
    #TODO: link here to all downloaded metadata. Rename resulting dataframes at will.
    SampleInfoDF = join(General.Basepath, 'SampleInfo.df')
    class bioGEOTRACES:
        metadataDF = join(General.Basepath, 'bioGEOTRACES_md.df')
        RawMetadataDir = join(General.Metadata,'GEOTRACES')
        CTDSensorTXT = join(RawMetadataDir, 'GEOTRACES_IDP2017_v1/ctd_sensor_data/ascii',
                            'GEOTRACES_IDP2017_v1_CTD_Sensor_Data.txt')
        DiscreteSampleTXT = join(RawMetadataDir, 'GEOTRACES_IDP2017_v2/discrete_sample_data/ascii',
                            'GEOTRACES_IDP2017_v2_Discrete_Sample_Data.txt')    
        SampleMeasurementsDF = join(General.Basepath, 'bioGEOTRACES_sample_meas.df')
    class ALOHA_BATS:
        metadataDF = join(General.Basepath, 'ALOHABATS_md.df')
        BATSRawMetadataDir = join(General.Metadata, 'BATS')
        BATSASCIIDir = join(BATSRawMetadataDir, 'ASCII')
        BATSBottle = join(BATSRawMetadataDir, 'bats_bottle.txt')
        BATSBottleHeader=59
        BATSSampleMeasurementsDF = join(General.Basepath, 'BATS_sample_meas.df')
        ALOHARawMetadataDir = join(General.Metadata, 'ALOHA')
        ALOHACTD = join(ALOHARawMetadataDir, 'ctd')
        ALOHAWater = join(ALOHARawMetadataDir, 'water')
        ALOHASampleMeasurementsDF = join(General.Basepath, 'ALOHA_sample_meas.df')
    class TARA:
        metadataDF = join(General.Basepath, 'TARA_md.df')
        RawMetadataDir = join(General.Metadata,'TARA')
        SeqContextXL = join(RawMetadataDir,'TARA_SAMPLES_CONTEXT_SEQUENCING_20170515.xlsx')
        Stations = join(RawMetadataDir, 'TARA_reg_stations.txt')
        CarbChemXL = join(RawMetadataDir, 'TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB_20170515.xlsx')
        NutrientsXL = join(RawMetadataDir, 'TARA_SAMPLES_CONTEXT_ENV-DEPTH-NUT_20170515.xlsx')
        HPLCXL = join(RawMetadataDir, 'TARA_SAMPLES_CONTEXT_ENV-DEPTH-HPLC_20170515.xlsx')
        XLHeaderLine = 17
        DepthSensorCSV = join(RawMetadataDir, 'TARA_ENV_DEPTH_SENSORS.tab')
        DepthSensorHeader = 2597
        SampleMeasurementsDF = join(General.Basepath, 'TARA_sample_meas.df')
    class United:
        metadataDF = join(General.Basepath, 'United_md.df')
        SampleMeasurementsDF = join(General.Basepath, 'United_sample_meas.df')
        
class RawFastq:
    # TODO: replace BasedDir and FastqDir in the following classes with where you place fastq
    # files for each of the datasets. Save sample sheet from ENA and have the SampleSheet
    # variable point to it.  
    class TARA:
        Basedir = join(General.Samples, 'TaraOceans')
        FastqDir = join(Basedir, 'PRJEB1787')
        SampleSheet = join(Basedir, 'PRJEB1787.txt')
    class bioGEOTRACES:
        Basedir = join(General.Samples, 'bioGEOTRACES')
        FastqDir = join(Basedir, 'PRJNA385854')
        SampleSheet = join(Basedir, 'PRJNA385854.txt')
    class ALOHA_BATS:
        Basedir = join(General.Samples, 'ALOHA-BATS')
        FastqDir = join(Basedir, 'PRJNA385855')
        SampleSheet = join(Basedir, 'PRJNA385855.txt')
        
class Mapping:
    class OM_RGC:
        IndexDir = join(General.Databases, 'OM-RGC')
        IndexFile = join(IndexDir,'OM-RGC')
        IndexFasta = join(IndexDir,'OM-RGC_seq.fasta')
        LengthsFile = join(Mapping.OM_RGC.IndexDir,'OM-RGC_seq.lens')
        MapDir = mkdirifnotexists(join(General.Tmppath,'OM-RGC_Mapping'))
        MapParams = dict(preset=MapPreset.SENSITIVE, report_alns=20, minins=0, maxins=500, 
                         no_mixed=False, no_discordant=False, dovetail=False, no_contain=False, 
                         no_overlap=False)
        ICRAParams = dict(max_mismatch=12, consider_lengths=True, epsilon=1e-6, \
                           max_iterations=30, min_bins=4, max_bins=100, min_reads=10, 
                           dense_region_coverage=60, length_minimum=300, \
                           length_maximum=2e5, use_theta=False, average_read_length=None, 
                           force_save_delta=True)
        ICRAUsage = OCEAN_GENES
        RemoveUnmapped = True 
        RemoveNotDelta=False
        DeltaThresh=0.999
        DeletePMP = True
        DeleteOldBAM=True
        
class Calling:
    class OM_RGC:
        CallDir = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_Call'))
        DbFasta = join(Mapping.OM_RGC.IndexDir,'OM-RGC_seq.fasta')
        FilterThreshold = 0.9
        mpileupParams = dict(min_base_qual = 15, max_depth_indel = 1000, 
                             max_depth = 100000, min_iReads = 2)
        
class SNP:
    class OM_RGC:
        InputDir = Calling.OM_RGC.CallDir
        GeneDFDir = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_GeneDFs'))
        OutDir = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_SNPAnalysis'))
        OutDirCollate = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_Collate'))
        OutDirCodons = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_Codons'))
        # Replace this with the path to which you ran eggNOG mapper on OM-RGC
        eggNOGMapper_outpath = join(Mapping.OM_RGC.IndexDir,'eggnog')
        GeneGroupCollateDBs = [join(eggNOGMapper_outpath,d) for d in \
                               ['KEGG_ko.dat','eggNOG_OGs.dat']]
        CacheDir = mkdirifnotexists(join(General.Tmppath, 'OM-RGC_Cache'))
        GeneLengths = join(Mapping.OM_RGC.IndexDir, 'OM-RGC_seq.lengths')
        OnlySNPs = True # Only SNPS or also indels
        QualThresh = 30 # Compound quality threshold
        MinSamples = 20 # Minimum number of samples
        MinVariants = 1 # Minimum number of variants per gene
        MinTotalVarSupport = 3 # Minimum number of variants to support a call
        MinMaf = 0.01 # Minimum minor allele frequency
        MinPosReads = 4 #Minimum coverage of reads per position to call SNP
        MinPercPoss = 60 #Minimal percent of positions complying with above number of reads
        MinGenes = 5 #Minimal number of genes per KEGG KO / eggNOG OG
    CodonProps = join(General.Basepath, 'Codon_costs.df')

        
class Analysis:
    class SNP:
        CountDir = join(General.Tmppath, 'OM-RGC_SNPCount')
        
#            
# class StrainGenes:
#     DBDir = join(General.Scratch, 'Data','Databases','GORG')
#     GenesTSV = join(DBDir, 'GORG_v1.tsv.gz')
#     SAGXLS = join(DBDir, 'gorg-tropics_sags_tableS2.xlsx')
#     MinStrainsPerGene = 20
#     MinStrainsPerGeneSubgrp = 10
#     AnalysisBase = mkdirifnotexists(join(General.Tmppath, 'GORG'))
#     FilteredSeqs = mkdirifnotexists(join(AnalysisBase, 'FiltSeqs'))
#     MSAs = mkdirifnotexists(join(AnalysisBase, 'MSAs'))
#     Trees = mkdirifnotexists(join(AnalysisBase, 'Trees'))
#     DFs = mkdirifnotexists(join(AnalysisBase, 'DFs'))
#     dNdSDir = mkdirifnotexists(join(AnalysisBase, 'dNdSDir'))