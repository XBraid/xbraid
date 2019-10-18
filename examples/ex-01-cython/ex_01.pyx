from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 

##
# See ex_01-setup.py for notes on installing and running
##

##
# Make braid_Vector and braid_App just PyObject* pointers
#   https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html
ctypedef PyObject* braid_Vector
ctypedef PyObject* braid_App


##
# Define your Python Braid Vector
cdef class PyBraid_Vector:
    '''
    Use a Python class to define your braid_Vector
    '''

    # This generic object can be set to a python object, like a numpy array
    cdef object value

    def __cinit__(self): 
        self.value = np.array((11.1,), dtype='float') 


    cdef void SetVectorPtr(self, braid_Vector *u_ptr):
        ''' 
        This member function assigns the pointer value u_ptr to point to
        this vector. 
        '''
        u_ptr[0] = <braid_Vector> self

##
# Define user App 
#
# Unlike with the C or Fortran Braid codes, the App structure is not as
# important, and does not need to be used. Instead, you can define "Global"
# Python Objects NOW that would normally go in App.  Then, you can use these
# objects below in any routine.  But, you may still want the App for "good style",
# so we include this example app object. 
#
##
cdef class PyBraid_App:
    '''
    Use a Python class to define your braid_App
    '''

    cdef int rank

    def __cinit__(self, _rank): 
        self.rank = _rank


##
# Import Cython Braid Wrappers 
# (only do after declaration of Braid Vector and App)
include "braid.pyx"


##
# Define user functions like Step, Init, ...
##

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    tstart = 0.0 # declare tstart and tstop as doubles
    tstop = 0.0
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    
    # Cast u as a PyBraid_Vector, and do time-stepping.  Note
    # how the value array is written in-place
    pyU = <PyBraid_Vector> u

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app

    #print("Step " + str(pyU.value))
    pyU.value[0] = 1./(1. + tstop-tstart)*pyU.value[0]
    #print("   Step Done" + str(pyU.value))

    return 0


cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    
    # Allocate new vector
    pyU = PyBraid_Vector() 
    
    if (t == 0.0):
        pyU.value[0] = 1.0
    else:
        pyU.value[0] = 0.456
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyU)
    
    # Set output, so that u_ptr points to pyU
    pyU.SetVectorPtr(u_ptr)

    #print("Init " + str(pyU.value))
    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    
    # Allocate new vector
    pyV = PyBraid_Vector()
    
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Assign pyU's value to clone pyV
    pyV.value[0] = pyU.value[0]
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyV)
    
    # Set output, so that u_ptr points to pyU
    pyV.SetVectorPtr(v_ptr)
    
    #print("Clone  " + str(pyU.value) + "  " + str(pyV.value))
    return 0
 
cdef int my_free(braid_App app, braid_Vector u):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Decrement the smart pointer
    Py_DECREF(pyU)

    del pyU 
    
    return 0

cdef int my_sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y):
    # Cast x and y as a PyBraid_Vector
    pyX = <PyBraid_Vector> x
    pyY = <PyBraid_Vector> y
    
    # Compute AXPY
    pyY.value[0] = alpha*pyX.value[0] + beta*pyY.value[0]

    #print("Sum " + str(pyY.value[0]) + ",  " + str(pyX.value[0]))
    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Compute norm 
    #  Note norm_ptr is a double array of size 1, and we index in at location [0]
    norm_ptr[0] = np.sqrt(np.dot(pyU.value, pyU.value))
    
    #print("Norm " + str(norm_ptr[0]) + ",  " + str(pyU.value[0]))
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app

    # Declare i as an integer, and the fill it in with the time index
    cdef int tindex
    braid_AccessStatusGetTIndex(status, &tindex)

    # Can get access to more information like
    # cdef double t
    # braid_AccessStatusGetT(status, &t)
    
    # Write this solution value to file
    filename = "%s.%04d.%03d"%("ex-01.out", tindex, pyApp.rank)
    f = open(filename, "w")
    f.write( "%.14e\n"%pyU.value);
    
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    
    #  Note size_ptr is an integer array of size 1, and we index in at location [0]
    size_ptr[0] = sizeof(double)
    
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u

    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
    
    # Pack buffer
    dbuffer[0] = pyU.value[0]

    braid_BufferStatusSetSize(status, sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    
    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
   
    # Allocate new vector
    pyU = PyBraid_Vector()
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyU)
    
    # Unpack information in dbuffer
    pyU.value[0] = dbuffer[0]
    
    # Set output, so that u_ptr points to pyU
    pyU.SetVectorPtr(u_ptr)
    
    #print("BufUnpack  " + str(pyU.value)) 
    return 0


def InitCoreApp():
    cdef braid_Core core
    cdef double tstart
    cdef double tstop
    cdef MPI.Comm comm = MPI.COMM_WORLD
    cdef int ntime
    cdef int rank
 
    ntime = 10
    tstart = 0.0
    tstop = tstart + ntime/2.
    rank = comm.Get_rank()

    # Declare app object
    pyApp = PyBraid_App(rank) 
    # Must increment smart pointer for app, or else Python will delete it 
    Py_INCREF(pyApp)

    # Initialize Braid
    braid_Init(comm.ob_mpi, comm.ob_mpi, tstart, tstop, ntime, 
               <braid_App> pyApp, my_step, my_init, my_clone, my_free, my_sum, 
               my_norm, my_access, my_bufsize, my_bufpack, 
               my_bufunpack, &core)

    # Store the Braid core inside of a Python-compatible object for return 
    pyCore = PyBraid_Core()
    pyCore.setCore(core)

    return pyCore, pyApp

def run_Braid(PyBraid_Core pyCore, PyBraid_App pyApp):
    
    # Set Braid options
    braid_SetMaxLevels(pyCore.getCore(), 2)
    braid_SetMaxIter(pyCore.getCore(), 10)

    # Run Braid
    braid_Drive(pyCore.getCore())
    
    # Destroy App, decrementing smart pointer for app
    Py_DECREF(pyApp)
    del pyApp 

    # Destroy Braid
    braid_Destroy(pyCore.getCore())


