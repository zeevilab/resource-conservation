from datetime import datetime
from gzip import open as opengz
import logging
import os
from os.path import expandvars, splitext, basename
from shutil import rmtree
from subprocess import check_output as sbcheck_output
import subprocess
import dill
from pandas import read_pickle
from pandas.io.pickle import to_pickle
log_ = logging.getLogger(__name__)
del logging
GZIP_EXT = '.gz'

def Write(outfile, datastruct, exceptiononmissingdir = False):
    dname, _ = os.path.split(outfile)
    if dname != '' and not os.path.exists(dname):
        assert not exceptiononmissingdir
        os.makedirs(dname)
    try:
        with open(outfile, 'wb') as pkl_file:
            dill.dump(datastruct, pkl_file, -1)
    except:
        to_pickle(datastruct, outfile)
        
def Load(infile, log = False):
    a = datetime.now()
    try:
        with open(infile, 'rb') as pkl_file: 
            res = dill.load(pkl_file)
    except:
        res = read_pickle(infile)
    b = datetime.now()
    if log:
        log_.info('Loaded '+ basename(infile) + ' in ' + str(b-a))
    return res
 
def mkdirifnotexists(pth):
    if not os.path.exists(pth):
        os.makedirs(pth)
    return pth
 
def chdirmkifnotexist(pth):
    if not os.path.exists(pth):
        os.makedirs(pth)
    os.chdir(pth)
    return pth

def split3way(fullpath):
    fdir, fname = os.path.split(fullpath)
    fname, ext = os.path.splitext(fname)
    return fdir, fname, ext

def shell_command(com, verbose = False, loglevel = None):
    com = expandvars(com)
    if loglevel is not None:
        if verbose:
            log_.log(loglevel, com)
        #CONSIDER:
        #    This doesn't support real-time logging (logging is
        #    done when the process finishes). See here for code 
        #    for the other option: 
        #    http://codereview.stackexchange.com/questions/6567/how-to-redirect-a-subprocesses-output-stdout-and-stderr-to-logging-module
        #    Though iter is probably better than while.
        errnout = ''
        try:
            errnout = sbcheck_output(com, shell = True, stderr = subprocess.STDOUT)
        except:
            if errnout != '':
                log_.log(loglevel, errnout)
            raise
        else:
            if errnout != '':
                log_.log(loglevel, errnout)
        return errnout
    else:
        log_.debug(com)
        return sbcheck_output(com, shell=True, stderr = subprocess.STDOUT)

def iswindows():
    return os.name == 'nt'

def rmrf(path, verbose = False):
    if iswindows():
        if verbose: print('Deleting ' + path)
        rmtree(path)
    else: #shutil has issues with some linux versions.
        shell_command('rm -rf ' + path, verbose)

def open_gz_indif(filepath):
    return opengz(filepath) if splitext(filepath)[1] == GZIP_EXT else open(filepath)

def rglob(basepath, wildcard):
    matches = []
    for root, _, filenames in os.walk(basepath):
        for filename in fnmatch.filter(filenames, wildcard):
            matches.append(os.path.join(root, filename))
    return matches

def tryrm(f):
    try:
        os.remove(f)
    except:pass
