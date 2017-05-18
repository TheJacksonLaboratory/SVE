from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules=[ Extension("fusion_utils",
              ["fusion_utils.pyx"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"])]

setup(
  name = 'fusion_utils',
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules,
  include_dirs=[np.get_include()]
)
