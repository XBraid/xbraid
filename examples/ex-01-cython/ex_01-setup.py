from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
# import numpy

##
# TODO for this branch
#
# Can we remove the BraidApp argument for the functions?  Maybe hold off on this.
# 
# Integrate these install comments automatically into the users manual.
# Perhaps have this be a separate directory inside of braid/examples/
#
# Design a few regression tests
#
# Can you make the install procedure more robust?
##

##
# This setup.py file has been tested on High Sierra with Homebrew, and Ubuntu LTS.
# This example requires Python 3 and the installation of Cython.
#
# To Install,
#
# 1) Make sure that library_dirs and include_dirs point to the 
#    location of "braid"
# 2) Type (using whatever install location you want)
#
#    $ python3 ex_01-setup.py install --prefix=$HOME/.local
#
#    Note that you may have to tweak the compilers and flags.
#    Some comments on this are below.
#
# To Run, 
#
# 1) Make sure that the install directory and the location of MPI4PY 
#    is in your PYTHONPATH, e.g.,
# 
#     export PYTHONPATH="$HOME/.local/lib/python3.6"
#     or perhaps, 
#     export PYTHONPATH="$HOME/.local/lib/python3.7/site-packages"
#
# 2) Type 
#    ipython 3
#    >>> import ex_01
#    >>> core = ex_01.InitCore()
#    >>> ex_01.run_Braid(core)
#
#    This should yield the same output as ex-01.c
#
# 3) For parallel runs, try
#
#    $$ mpirun -np K  python3 ex_01_run.py
## 

##
# Other notes:
#  1) Some systems may need to find Numpy headers, which are located in 
#     include_dirs=["../braid", numpy.get_include()],
#  2) Some compilers may require "-fPIC" to be added to extra_compile_args
#
##


os.environ["CC"] = "mpicc"
os.environ["LDSHARED"] = "mpicc -shared"    # Comment out for High Sierra with Homebrew

ex_01_extension = Extension(
    name="ex_01",
    sources=["ex_01.pyx"],
    libraries=["braid"],
    library_dirs=["../../braid"],
    include_dirs=["../../braid"], 
    extra_compile_args=["-Wno-incompatible-pointer-types", "-Wno-unused-function"] 
)
setup(
    name="ex_01",
    ext_modules=cythonize([ex_01_extension], language_level = "3")
)

