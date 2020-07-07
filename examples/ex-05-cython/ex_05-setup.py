from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import sys
# import numpy

braid_dir = "../../braid"

##
# This setup.py file has been tested on High Sierra with Homebrew, and Ubuntu LTS.
# This example requires Python 3 and the installation of Cython.
#
# To Install,
#
# 1) Make sure that braid_dir (defined above) points to the location of "braid"
#    This is used by library_dirs and include_dirs below
#
# 2) Type (using whatever install location you want)
#
#    $ python3 ex_05-setup.py install --prefix=$HOME/.local
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
#    python3
#    >>> import ex_05
#    >>> core, app = ex_05.InitCoreApp()
#    >>> ex_05.run_Braid(core, app)
#
#    Output:
#
#
## 3) For parallel runs, try
#
#    $$ mpirun -np K  python3 ex_05_run.py
## 

##
# Other notes:
#  1) Some systems may need to find Numpy headers, which are located in 
#     include_dirs=["../braid", numpy.get_include()],
#  2) Some compilers may require "-fPIC" to be added to extra_compile_args
#
##


os.environ["CC"] = "mpicc"
if sys.platform != 'darwin':
    os.environ["LDSHARED"] = "mpicc -shared"    # Needed by Ubuntu LTS, leave out for Mac

ex_05_extension = Extension(
    name="ex_05",
    sources=["ex_05.pyx"],
    libraries=["braid"],
    library_dirs=[braid_dir],
    include_dirs=[braid_dir], 
    extra_compile_args=["-Wno-incompatible-pointer-types", "-Wno-unused-function" ] 
)
setup(
    name="ex_05",
    ext_modules=cythonize([ex_05_extension], language_level = "3")
)

