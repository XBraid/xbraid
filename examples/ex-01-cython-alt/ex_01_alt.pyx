from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from math import sqrt
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 


##
# Define your Python Braid Vector as a C-struct
#
# Although this is a scalar problem, we still make "v" a <double *> so that
# this example shows how to create, allocats, and deallocates C-arrays 
cdef struct _braid_Vector_struct:
    double* v
ctypedef _braid_Vector_struct my_Vector
ctypedef _braid_Vector_struct *braid_Vector


##
# Define your App as a C-struct
cdef struct _braid_App_struct:
    int rank
ctypedef _braid_App_struct my_App
ctypedef _braid_App_struct *braid_App


##
# Import Cython Braid Wrappers 
# (only do after declaration of Braid Vector and App)
include "braid.pyx"


##
# Define user functions like Step, Init, ...
##

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    cdef double tstart
    cdef double tstop
    tstart = 0.0
    tstop = 5.0
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    
    # Now, we show how to make u.v look like a NumPy array
    cdef double[:] v_view = <double[:1]> u.v    # v_view is a memory view of length-1 array u.v
    v_np = np.asarray(v_view)                   # we can cast memory views as numpy arrays

    # Note this operation is done with a Numpy array object, but writes to the memory in u.v
    v_np[0] = 1./(1. + tstop-tstart)*v_np[0]

    return 0

cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    
    # Manually allocate a new vector u
    cdef my_Vector* u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    u.v = <double*>PyMem_Malloc(1*sizeof(double))
    
    # Create Numpy wrapper around u.v
    cdef double[:] v_view = <double[:1]> u.v
    v_np = np.asarray(v_view)

    if (t == 0.0):
        v_np[0] = 1.0
    else:
        v_np[0] = 0.456

    u_ptr[0] = u        # Weird Cython way of de-referencing u_ptr
    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    
    # Manually allocate a new vector v
    cdef my_Vector* v
    v = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    v.v = <double*>PyMem_Malloc(1*sizeof(double))
    
    # Create Numpy wrapper around v.v
    cdef double[:] vv_view = <double[:1]> v.v
    v_np = np.asarray(vv_view)

    # Create Numpy wrapper around u.v
    cdef double[:] uv_view = <double[:1]> u.v
    u_np = np.asarray(uv_view)

    v_np[:] = u_np[:]

    v_ptr[0] = v        # Weird Cython way of de-referencing v_ptr
    return 0

cdef int my_free(braid_App app, braid_Vector u):
    PyMem_Free(u.v)
    PyMem_Free(u)
    return 0

cdef int my_sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y):
    # Create Numpy wrapper around x.v
    cdef double[:] xv_view = <double[:1]> x.v
    x_np = np.asarray(xv_view)

    # Create Numpy wrapper around y.v
    cdef double[:] yv_view = <double[:1]> y.v
    y_np = np.asarray(yv_view)

    # Must be carefuly to write y_np IN-PLACE  with [:] notation !!
    y_np[:] = alpha*x_np + beta*y_np
    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    
    # Create Numpy wrapper around u.v
    cdef double[:] uv_view = <double[:1]> u.v
    u_np = np.asarray(uv_view)
    
    norm_ptr[0] = sqrt(np.dot(u_np, u_np))  # Weird Cython way of de-referencing norm_ptr
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    cdef double t
    
    # Create Numpy wrapper around u.v
    cdef double[:] uv_view = <double[:1]> u.v
    u_np = np.asarray(uv_view)

    braid_AccessStatusGetT(status, &t)
    print(u_np[0], t)
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    
    size_ptr[0] = sizeof(double) # Weird Cython way of de-referencing size_ptr
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer
    
    # Create Numpy wrapper around u.v
    cdef double[:] uv_view = <double[:1]> u.v
    u_np = np.asarray(uv_view)

    dbuffer[0] = u_np[0]

    braid_BufferStatusSetSize(status, sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):

    cdef double *dbuffer = <double*> buffer

    # Manually allocate a new vector u
    cdef my_Vector* u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    u.v = <double*>PyMem_Malloc(1*sizeof(double))
    
    # Create Numpy wrapper around u.v
    cdef double[:] v_view = <double[:1]> u.v
    v_np = np.asarray(v_view)

    # Create Numpy wrapper around dbuffer
    cdef double[:] buf_view = <double[:1]> dbuffer
    buf_np = np.asarray(buf_view)

    v_np[:] = buf_view[:]

    u_ptr[0] = u # Weird Cython way of de-referencing u_ptr
    return 0


cdef class PyBraid_App:
    '''
    This class is a simple wrapper about the C-struct braid_App
    This class allows us to return a braid_App C-struct back into Python,
    
    Note, that you cannot pass C-structs like braid_App directly back to
    Python
    '''
    cdef braid_App app

    def __init__(self): 
        pass

    cdef setApp(self, braid_App _app):
        self.app = _app
    
    cdef braid_App getApp(self):
        return self.app


##
# Sets up Braid for running, initialize the App and Core
def InitCoreApp():

    cdef braid_Core core
    cdef my_App *app            # Declare App as a C-struct
    cdef double tstart
    cdef double tstop
    cdef MPI.Comm comm = MPI.COMM_WORLD
    cdef int ntime
    cdef int rank
 
    ntime = 10
    tstart = 0.0
    tstop = tstart + ntime/2.
    rank = comm.Get_rank()
    app = <my_App*>PyMem_Malloc(sizeof(my_App))
    app.rank = rank

    braid_Init(comm.ob_mpi, comm.ob_mpi, tstart, tstop, ntime, app, my_step,
            my_init, my_clone, my_free, my_sum, my_norm, my_access, my_bufsize,
            my_bufpack, my_bufunpack, &core)
    
    # Store the Braid core inside of a Python-compatible object for return 
    pyCore = PyBraid_Core()
    pyCore.setCore(core)

    # We have to do a similar trick for the app, 
    pyApp = PyBraid_App()
    pyApp.setApp(app)

    return pyCore, pyApp


##
# Run Braid with App and Core
def run_Braid(PyBraid_Core pyCore, PyBraid_App pyApp):
    
    # Set Braid options
    braid_SetMaxLevels(pyCore.getCore(), 2)
    braid_SetMaxIter(pyCore.getCore(), 10)

    # Run Braid
    braid_Drive(pyCore.getCore())
    
    # Destroy App C-Struct
    PyMem_Free(pyApp.getApp())

    # Destroy Braid Core C-Struct
    braid_Destroy(pyCore.getCore())


