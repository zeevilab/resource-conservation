from os.path import join, exists
from Bio import SeqIO
from shutil import copy
from glob import glob
from time import sleep
import os
from pandas.core.reshape.concat import concat
from pandas.io.parsers import read_csv
from _collections import defaultdict
from pandas.io.pickle import read_pickle
from lib.Utils import split3way, rmrf, chdirmkifnotexist, shell_command, Write
from Data.config import Annotate

def translate_db(db_in, db_out):
    with open(db_out, 'w') as fout:
        for rec in SeqIO.parse(db_in, 'fasta'):
            newrec = rec
            newrec.seq = rec.seq.translate(table=11)[:-1]
            newrec.description = ''
            SeqIO.write(newrec, fout, 'fasta')

def split_database(db_in, dirout, chunk_size):
    with open(db_in) as fin:
        chunknum = 0
        fout = open(join(dirout, '{}_{:04d}.fa'.format(split3way(db_in)[1],chunknum)),'w')
        culchunk = 0
        for rec in SeqIO.parse(fin, 'fasta'):
            culchunk += len(rec)
            if culchunk > chunk_size:
                fout.close()
                chunknum += 1
                culchunk = 0
                fout = open(join(dirout, '{}_{:04d}.fa'.format(split3way(db_in)[1],chunknum)),'w')
            SeqIO.write(rec, fout, 'fasta')
        fout.close()

def eggnog_map_one(f, cpu=8):
    for db_f in ['eggnog.db','eggnog_proteins.dmnd']:
        if not exists(join('/dev/shm', db_f)):
            copy(join(Annotate.OM_RGC.EmapperDir, 'data', db_f), join('/dev/shm', db_f))
        elif os.stat(join('/dev/shm', db_f)).st_size != os.stat(join(Annotate.OM_RGC.EmapperDir, 
                                                                     'data', db_f)).st_size:
            sleep(200) #wait for file to copy
    tmpdir = join(Annotate.OM_RGC.AnnotDir, 'tmp', split3way(f)[1])
    rmrf(tmpdir)
    chdirmkifnotexist(tmpdir)
    shell_command("{} {} -i {} --output {} -m diamond --cpu {} --override --data_dir /dev/shm --temp_dir {} --no_file_comment"\
                  .format(Annotate.OM_RGC.Py27, Annotate.OM_RGC.EmapperPy, f, 
                          join(Annotate.OM_RGC.AnnotDir, split3way(f)[1]), cpu, 
                          tmpdir), verbose=True)

def parseeggnog():
    df_f = join(Annotate.OM_RGC.EggnogDir, 'OM-RGC_annotations.df')
    if exists(df_f):
        df = read_pickle(df_f)
    else:
        df = concat([read_csv(f, sep='\t', header=None) \
                     for f in sorted(glob(join(Annotate.OM_RGC.AnnotDir, '*.annotations')))])
        df.columns = ['GeneID', 'seed_eggNOG_ortholog', 'seed_ortholog_evalue',
                      'seed_ortholog_score', 'Predicted_taxonomic_group', 'Predicted_protein_name',
                      'GeneOntology', 'EC_number', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
                      'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction',
                      'tax_scope', 'eggNOG_OGs', 'bestOG', 'COG_Functional_Category',
                      'eggNOG_description']
        df = df.set_index('GeneID')
        df.to_pickle(df_f)
    
    for col in df.columns[3:19]:
        if col == 'bestOG':
            continue
        if exists(join(Annotate.OM_RGC.EggnogDir, col + '.dat')):
            continue
        ret = defaultdict(list)
        locser = df[col].dropna()
        for nm in locser.index:
            if col == 'eggNOG_OGs':
                terms = set([x.split('@')[0] for x in locser[nm].split(',')])
            elif col == 'COG_Functional_Category':
                terms = set([x for x in locser[nm]])
            else:
                terms = set(locser[nm].split(','))
            for term in terms:
                ret[term].append(nm)
        Write(join(Annotate.OM_RGC.EggnogDir, col + '.dat'), ret)
        print(col)
    