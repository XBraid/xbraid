from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
# import numpy

#
# See ex-01-cython/ex_01-setup.py for installation instructions
#

os.environ["CC"] = "mpicc"
#os.environ["LDSHARED"] = "mpicc -shared"    # Comment out for High Sierra with Homebrew

ex_01_alt_extension = Extension(
    name="ex_01_alt",
    sources=["ex_01_alt.pyx"],
    libraries=["braid"],
    library_dirs=["../../braid"],
    include_dirs=["../../braid"],
    extra_compile_args=["-Wno-incompatible-pointer-types", "-Wno-unused-function"] 
)
setup(
    name="ex_01_alt",
    ext_modules=cythonize([ex_01_alt_extension], language_level = "3")
)

