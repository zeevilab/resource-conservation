from Data.Metadata.TARA_metadata import getbasemetadata, get_measurements_TARA
from Data.config import Biodata, RawFastq, Annotate
from Data.Metadata.bioGEOTRACES_metadata import extract, get_measurements_GEOTRACES,\
    get_measurements_BATS, get_measurements_HOT
from Data.Metadata.KEGGAnnotate import translate_db, split_database, eggnog_map_one, parseeggnog
from glob import glob
from os.path import join, exists
from lib.Utils import split3way
from Data.Metadata.UniteMetadata import unite_sampledata, unite_measurements

def metadata():
    # Get Tara metadata
    getbasemetadata()
    get_measurements_TARA()
    # Get bioGEOTRACES metadata
    extract(RawFastq.bioGEOTRACES, Biodata.bioGEOTRACES)
    get_measurements_GEOTRACES(depth_tolerance=5)
    # Get HOT/BATS metadata
    extract(RawFastq.bioGEOTRACES, Biodata.bioGEOTRACES)
    get_measurements_BATS()
    get_measurements_HOT()
    # Unite metadata from all sources
    unite_sampledata()
    unite_measurements()
    
def annotate():
    translate_db(Annotate.OM_RGC.IndexFasta, Annotate.OM_RGC.IndexFaa)
    # Split the OM-RGC database to chuncks to reduce runtime
    split_database(Annotate.OM_RGC.IndexFaa, Annotate.OM_RGC.SplitDir, 5e7)
    # Run eggnog mapper on each one of the chunks
    # TODO: IMPORTANT! Wrap this for loop with your hpc job submission pipeline 
    for fname in glob(join(Annotate.OM_RGC.SplitDir, '*.fa')):
        if exists(join(Annotate.OM_RGC.AnnotDir, split3way(fname)[1]) + '.emapper.annotations'):
            continue
        eggnog_map_one(fname, 8)
    # Parse the resulting annotation database and create a dictionary file for each annotation
    # This is then used in the collate stage (SNP pipeline)
    parseeggnog()
    
def main():
    metadata()
    annotate()

if __name__ == '__main__':
    main()