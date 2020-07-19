import pickle
# import csv
from datetime import datetime
# import errno
# import fnmatch
# from glob import glob
from gzip import open as opengz
# from itertools import zip_longest
import logging
import os
from os.path import expandvars, splitext, join, basename
# import platform
# from random import choice
# import re
from shutil import rmtree
# import string
from subprocess import check_output as sbcheck_output
import subprocess
# import sys
# from timeit import default_timer as timer
# 
import dill
from pandas import read_pickle
# from pandas.core.frame import DataFrame
from pandas.io.pickle import to_pickle
# 
# from addloglevels import addloglevels, LOW_DEBUG
# from functools import wraps
# from logging import DEBUG
# 
# try:
#     import xlwt
#     import xlrd
# except:
#     print ('You need to download and install xlwt, xlrd')
# 
# CAT = 'cat'
# ZCAT = 'zcat'
# 
log_ = logging.getLogger(__name__)
del logging
# addloglevels()
GZIP_EXT = '.gz'
# 
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
# 
# def write_iterable(outfile, liststruct, exceptiononmissingdir = False):
#     dname, _ = os.path.split(outfile)
#     if dname != '' and not os.path.exists(dname):
#         assert not exceptiononmissingdir
#         os.makedirs(dname)
#     with open(outfile, 'wb') as pkl_file: 
#         for item in liststruct:
#             dill.dump(item, pkl_file, -1)
# 
# class buffered_writer:
#     def __init__(self,outdesc, buffsize):
#         self.outdesc = outdesc
#         self.buffsize = buffsize
#         self.buffer = []
#         
#     def write(self, strct):
#         self.buffer.append(strct)
#         if len(self.buffer) == self.buffsize:
#             self.flush()
#     
#     def flush(self):
#         dill.dump(self.buffer, self.outdesc, -1)
# 
# def load_iterable(infile):
#     try:
#         with open(infile, 'rb') as pkl_file:
#             while True:
#                 yield dill.load(pkl_file)
#     except EOFError:
#         raise StopIteration
#     
# def load_iterable_fromdesc(desc):
#     try:
#         while True:
#             yield pickle.load(desc)
#     except EOFError:
#         raise StopIteration
# 
# class TimerWrapper():
#     def __init__(self):
#         self.elapsed_sec = 0
#         self.is_running = False
#         self._curstart = None
#         
#     def start(self):
#         assert not self.is_running
#         self.is_running = True
#         self._curstart = timer()
#         
#     def stop(self):
#         stoppedat = timer()
#         assert self.is_running
#         self.is_running = False
#         self.elapsed_sec += (stoppedat - self._curstart)
#         self._curstart = None
#         
#     def reset(self):
#         self.elapsed_sec = 0
#         self.is_running = False
#         self._curstart = None
# 
# def none_iter():
#     while True:
#         yield None
# 
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
 
# def LoadCSV(infile, delim = ','):
#     with open(infile, 'rb') as instrm:
#         lol = list(csv.reader(instrm, delimiter = delim))
#     return lol
# 
# def LoadCSVIter(infile, delim = ','):
#     with open(infile, 'rb') as instrm:
#         rdr = csv.reader(instrm, delimiter = delim)
#         for line in rdr:
#             yield line
#             
# def LoadXLS(infile):
#     results = {}
#     workbook = xlrd.open_workbook(infile)
#     for ws_name in workbook.sheet_names():
#         worksheet = workbook.sheet_by_name(ws_name)
#         cur = []
#         for r in range(worksheet.nrows):
#             cur.append(worksheet.row(r))
#         results[ws_name] = cur
#     return results
#             
# def SaveCSV(outfile, rows, delim = ','):
#     with open(outfile, 'wb') as outstrm:
#         csv.writer(outstrm, delimiter =  delim).writerows(rows)
#         
def mkdirifnotexists(pth):
    if not os.path.exists(pth):
        os.makedirs(pth)
    return pth
 
def chdirmkifnotexist(pth):
    if not os.path.exists(pth):
        os.makedirs(pth)
    os.chdir(pth)
    return pth
# 
# def floatifpos(val):
#     try:
#         return float(val)
#     except:
#         return val
#     
# def outlog(string, fh = sys.stderr):
#     fh.write("{} {}\n".format(str(datetime.now())[:-7], string))
# 
# def unitecsvs(path, extension = 'csv', delimiter = ',', name = ''):
#     finalname = join(path,name+"compiled.xls")
#     if os.path.exists(finalname):
#         os.remove(finalname)
#     wb = xlwt.Workbook()
#     for filename in sorted(glob(join(path,"*.%s" % extension))):
#         (_, f_name) = os.path.split(filename)
#         (f_short_name, _) = os.path.splitext(f_name)
#         f_short_name = f_short_name[-30:]
#         ws = wb.add_sheet(f_short_name)
#         with open(filename, 'rb') as fle:
#             spamReader = csv.reader(fle, delimiter = delimiter)
#             for rowx, row in enumerate(spamReader):
#                 for colx, value in enumerate(row):
#                     ws.write(rowx, colx, floatifpos(value))
#         os.remove(filename)
#     wb.save(finalname)
# 
# #taken from http://stackoverflow.com/questions/5251330/python-function-to-make-arbitrary-strings-valid-filenames
# def removeDisallowedFilenameChars(filename):
#     return re.sub(r"[\/\*\?\"\<\>\|]", '', filename)
# 
def split3way(fullpath):
    fdir, fname = os.path.split(fullpath)
    fname, ext = os.path.splitext(fname)
    return fdir, fname, ext
# 
# def rmifexists(p):
#     #OMG, this is so pythonic. 
#     #http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
#     try:
#         os.remove(p)
#     except OSError as e: # this would be "except OSError, e:" before Python 2.6
#         if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
#             raise # re-raise exception if a different error occured)
#         
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
# 
# def wcl(fpath):
#     if iswindows():
#         return sum(1 for _ in open(fpath, 'rb'))
#     else:
#         return int(shell_command('wc -l ' + fpath).split()[0])
# 
def iswindows():
    return os.name == 'nt'
# 
# def ismac():
#     return os.name == 'posix' and platform.system() == 'Darwin'
# 
# def is_es_comp():
#     return os.path.exists('/Users/eran')
# 
# 
# def is_is_comp():
#     return os.path.exists('/home/ira/')
# 
# 
# def is_ubuntu():
#     return os.name == 'posix' and platform.dist()[0] == 'Ubuntu'
# 
def rmrf(path, verbose = False):
    if iswindows():
        if verbose: print('Deleting ' + path)
        rmtree(path)
    else: #shutil has issues with some linux versions.
        shell_command('rm -rf ' + path, verbose)
# 
# #http://stackoverflow.com/questions/1835018/python-check-if-an-object-is-a-list-or-tuple-but-not-string/1835259#1835259    
# def is_sequence(arg):
#     return (not hasattr(arg, "strip") and
#             hasattr(arg, "__getitem__") or
#             hasattr(arg, "__iter__"))
#     
# def chmod_r(path, mode):
#     for root, dirs, files in os.walk(path):
#         for d in dirs:
#             try: os.chmod(os.path.join(root, d), mode)
#             except: pass
#         for f in files: 
#             try: os.chmod(os.path.join(root, f), mode)
#             except: pass
#             
def open_gz_indif(filepath):
    return opengz(filepath) if splitext(filepath)[1] == GZIP_EXT else open(filepath)
# 
# def seqioparse_gz_indif(filepath):
#     #BioPython is a problematic install, and we don't want it to contaminate Utils.
#     from Bio import SeqIO 
#     return SeqIO.parse(opengz(filepath), 'fastq') if splitext(filepath)[1] == GZIP_EXT else \
#         SeqIO.parse(open(filepath), 'fastq')
#         
# def md5checksum(filepath):
#     return shell_command('md5sum ' + filepath).split()[0]
# 
# def rglob(basepath, wildcard):
#     matches = []
#     for root, _, filenames in os.walk(basepath):
#         for filename in fnmatch.filter(filenames, wildcard):
#             matches.append(os.path.join(root, filename))
#     return matches
# 
# #http://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression
# def merge_dicts(*dict_args):
#     '''
#     Given any number of dicts, shallow copy and merge into a new dict,
#     precedence goes to key value pairs in latter dicts.
#     '''
#     result = {}
#     for dictionary in dict_args:
#         result.update(dictionary)
#     return result
# 
# def name_value_pair_str_2_dict (s):
#     res = {}
#     for item in s:
#         name, var = item.partition("=")[::2]
#         res[name] = var
#     return re
# 
# tax = ('k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
# taxheaders = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
# 
# def isint(a):
#     try:
#         _ = int(a)
#         return True
#     except:
#         return False
#     
# #https://docs.python.org/2/library/itertools.html#recipes
# def grouper(iterable, n, fillvalue=None):
#     "Collect data into fixed-length chunks or blocks"
#     # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
#     args = [iter(iterable)] * n
#     return zip_longest(fillvalue=fillvalue, *args)
# 
# def dataFrameSplitRowsByFuncOnColumn(df,columntosplit,func=list):
#     """
# *****Example1*****
# a = pandas.DataFrame([[1,4,'A'],[2,5,'$$'],[3,6,'$G']],columns=['one','two','three'])
#     one  two three
# 0    1    4     A
# 1    2    5    $$
# 2    3    6    $G
# 
# dataFrameSplitRowsByFuncOnColumn(a,'three')
#    one  two three
# 0    1    4     A
# 1    2    5     $
# 2    2    5     $
# 3    3    6     $
# 4    3    6     G
# 
# *****Example2*****
# a = pandas.DataFrame([[1,4,'A'],[2,5,'$,$'],[3,6,'$,G']],columns=['one','two','three'])
#    one  two three
# 0    1    4     A
# 1    2    5   $,$
# 2    3    6   $,G
# 
# dataFrameSplitRowsByFuncOnColumn(a,'three',func=lambda x: x.split(","))
#    one  two three
# 0    1    4     A
# 1    2    5     $
# 2    2    5     $
# 3    3    6     $
# 4    3    6     G
#     """
#     columnnames=df.columns.values.tolist()
#     columnidx=columnnames.index(columntosplit)
#     columnnames.remove(columntosplit)
#     iterdf = df.set_index(columnnames)
#     res = DataFrame(iterdf[columntosplit].apply(func).tolist(),index=iterdf.index).stack()
#     columnnames.insert(columnidx,0)
#     res = res.reset_index()[columnnames].rename(columns={0:columntosplit})
#     return res 
# 
# import inspect
# def in_test():
#     current_stack = inspect.stack()
#     for stack_frame in current_stack:
#         if stack_frame is None:
#             continue
#         if stack_frame[4] is None:
#             continue
#         for program_line in stack_frame[4]: 
#             if "unittest" in program_line:      
#                 return True
#         if "/tests/test_" in stack_frame[1]:
#                 return True
#     return False
# 
# def calcFilenameForCacheOnDisk(decorator_self, **kwargs):
#     if decorator_self.basePath is not None:
#         path = os.path.join(decorator_self.basePath, decorator_self.func_name)
#         mkdirifnotexists(path)
#     else:
#         path = ''
#     filename=os.path.join(path,decorator_self.filename % kwargs)
#     return filename
# 
# def run_criterion_for_CacheOnDisk(decorator_self, **kwargs):
#     filename = calcFilenameForCacheOnDisk(decorator_self, **kwargs)
#     shouldRunFunction = not os.path.exists(filename)
#     if  decorator_self.force:
#         shouldRunFunction = True
#     if 'mustBeCached_ifNotThen' in kwargs:
#         shouldRunFunction = False
#     return shouldRunFunction  
# 
# class cacheOnDisk(object):
#     def __init__(self,filename='self.dat',basePath=None,force=False, redoOnNone=False, alsoInDebug=False):
#         self.basePath = basePath
#         self.filename = filename
#         self.force = force
#         self.inUnitTest = in_test()
#         self.redoOnNone = redoOnNone
#         self.alsoInDebug = alsoInDebug
# 
#     def __call__(self, function):
#         decorator_self = self
#         decorator_self.func_name = function.func_name
#         def wrappee( *args, **kwargs):
#             if not decorator_self.inUnitTest or self.alsoInDebug:
#                 filename = calcFilenameForCacheOnDisk(decorator_self, **kwargs)
#                 if not run_criterion_for_CacheOnDisk(decorator_self, **kwargs): 
#                     print ("Cache on disk loading %s" % filename)
#                     res = Load(filename)
#                     print ("Done.")
#                     if res is not None or not decorator_self.redoOnNone:
#                         return res
#                 if 'mustBeCached_ifNotThen' in kwargs:
#                     exec **kwargs['mustBeCached_ifNotThen'] in globals(), locals()
#                     
#             res=function(*args, **kwargs)
#             if not decorator_self.inUnitTest:
#                 print ("Cache on disk saving %s" % filename)
#                 Write(filename,res)
#                 print ("Done.")
#             return res
#         wrappee.original_func_name = function.func_name
#         wrappee.original_module_name = function.__module__
#         wrappee.original_decorator_self = decorator_self
#         return wrappee
#     
# class dontRunIfFileExists(object):
#     def __init__(self,filename,force=False):
#         self.filename = filename
#         self.force = force
#         self.inUnitTest = in_test()
#     def __call__(self, function):
#         decorator_self = self
#         def wrappee( *args, **kwargs):
#             try:
#                 filename=decorator_self.filename % kwargs
#             except:
#                 filename=decorator_self.filename % args
#             if os.path.exists(filename) and not decorator_self.force and not decorator_self.inUnitTest:
#                 return None
#             return function(*args, **kwargs)           
#         return wrappee
# 
# def catcom(do_zcat, infilelist, outfile,append=False, loglevel = LOW_DEBUG):
#     com = '{} {} > {}'.format(ZCAT if do_zcat else CAT, ' '.join(sorted(infilelist)), \
#                               outfile)
#     if append:
#         if os.path.exists(outfile):
#             com = '{} {} >> {}'.format(ZCAT if do_zcat else CAT, ' '.join(sorted(infilelist)), outfile)  
#     
#     shell_command(com, loglevel = loglevel)
# ops = None
# def randstr(n = 12):
#     global ops
#     if ops is None:
#         ops = string.ascii_letters + string.digits
#     return ''.join(choice(ops) for _ in range(n))
# 
# def ftplistdir(ftp_url, login=None, password=None):
#     import ftplib
#     from urllib.parse import urlparse
#     ftp = ftplib.FTP(urlparse(ftp_url).netloc)
#     if login is None:
#         ftp.login()
#     else:
#         ftp.login(login, password)
#     ftp.cwd(urlparse(ftp_url).path)
#     try:
#         files = ftp.nlst()
#     except ftplib.error_perm as resp:
#         if str(resp) == "550 No files found":
#             return []
#         else:
#             raise
#     return files
# 
# def timeit(log_, text = None, level = DEBUG):
#     def _inner(method):
#         '''A decorator which logs the time it took a method to run into debug'''
#         @wraps(method)
#         def timed(*args, **kw):
#             ts = datetime.now()
#             result = method(*args, **kw)
#             te = datetime.now()
#             log_.log(level, ('Method {} took {}'.format(method.__name__, te-ts) \
#                              if text is None else text.format(method.__name__, te-ts)))
#             return result
#         return timed
#     return _inner
# 
def tryrm(f):
    try:
        os.remove(f)
    except:pass
