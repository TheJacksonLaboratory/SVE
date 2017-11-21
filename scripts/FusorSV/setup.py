#!/usr/env/bin/python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import pysam
import numpy
import glob
import os

cythonize('fusion_utils.pyx')
extensions = [Extension('fusion_utils',
                        sources=['fusion_utils.pyx'],
                        include_dirs=pysam.get_include()+[numpy.get_include()],
                        define_macros=pysam.get_defines(),
                        extra_compile_args=['-ffast-math'])]

setup(
    name = 'fusorsv',
    author='Wan-Ping Lee',
    author_email='wan-ping.lee@jax.org',
    license='GPL 3 License',
    description='SV calling data fusion framework',
    classifiers=['Intended Audience :: Developers',
                 'License :: GPL 3 License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Cython',
                 'Programming Language :: C',
                 'Operating System :: POSIX',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    cmdclass = { 'build_ext': build_ext },
    ext_modules = extensions,
    scripts    = ['FusorSV.py'])
