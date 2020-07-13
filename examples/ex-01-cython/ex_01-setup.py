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
#    python3
#    >>> import ex_01
#    >>> core, app = ex_01.InitCoreApp()
#    >>> ex_01.run_Braid(core, app)
#
#    Print output with
#    >>> import os; os.system("cat ex-01.out.00*")
#    1.00000000000000e+00
#    6.66666666666667e-01
#    4.44444444444444e-01
#    2.96296296296296e-01
#    1.97530864197531e-01
#    1.31687242798354e-01
#    8.77914951989026e-02
#    5.85276634659351e-02
#    3.90184423106234e-02
#    2.60122948737489e-02
#    1.73415299158326e-02
#    0
#
## 3) For parallel runs, try
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
if sys.platform != 'darwin':
    os.environ["LDSHARED"] = "mpicc -shared"    # Needed by Ubuntu LTS, leave out for Mac

ex_01_extension = Extension(
    name="ex_01",
    sources=["ex_01.pyx"],
    libraries=["braid"],
    library_dirs=[braid_dir],
    include_dirs=[braid_dir], 
    extra_compile_args=["-Wno-incompatible-pointer-types", "-Wno-unused-function"] 
)
setup(
    name="ex_01",
    ext_modules=cythonize([ex_01_extension], language_level = "3")
)

