import os, sys
import numpy
from numpy.distutils.core import setup, Extension

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

#"mmio.c", "mm.cc", "mshift.cc", "eigen_numpy.cc", 
module = Extension('flc', 
                   sources = ["flc.cc", "mmio.c", "eigen_numpy.cc"],
                   library_dirs = ["./libs/_install/lib"],
                   libraries = ["boost_numpy", "boost_python", "glog"],
                   extra_compile_args = ["-w", "-std=gnu++11"],
                   include_dirs = [numpy.get_include(), "/usr/include/eigen3", "./libs/_install/include"] )

setup(name = 'flc',
      version = '0.1',
      ext_modules=[module])
