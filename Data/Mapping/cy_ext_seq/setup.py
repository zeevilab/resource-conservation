from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy as np
ext_modules = [Extension("gem_qa", ["gem_qa.pyx"], include_path = [np.get_include(), './'], include_dirs = [np.get_include(), './'], 
                         language = 'c++', extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]),\
               Extension("CSeqDict", ["CSeqDict.pyx"], include_path = ['./'], include_dirs = [], language = 'c++', extra_compile_args=["-std=c++11"],
                         extra_link_args=["-std=c++11"]),
               Extension("perc_id", ["perc_id.pyx"], include_path = [np.get_include(), './'], include_dirs = [np.get_include(), './'], \
                         language = 'c++', extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]),
               Extension("sam2pmp_helper", ["sam2pmp_helper.pyx"], include_path = [np.get_include(), './'], include_dirs = [np.get_include(), './'], 
                         language = 'c++', extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"])]
setup(cmdclass={'build_ext': build_ext}, ext_modules=ext_modules)
