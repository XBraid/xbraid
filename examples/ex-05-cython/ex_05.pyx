from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 
from scipy.sparse import spdiags, eye, linalg, csr_matrix
from sys import exit
from invert_sparse_mat_splu import invert_sparse_mat_splu 

#
# Example:       ex_05.pyx
#
# Interface:     Cython
# 
# Requires:      Python 3, Cython, C-language support     
#
# Description:   Solve the 1D advection-diffusion equation 
#                   u_t = u_x + eps*u_xx 
#                
#                on x in [0, 1], with t in [0, tstop=1.0]
#
#                Initial condition is a step function.
#
#                Boundary conditions are 
#                u(x=0, t) = 1
#                u(x=1, t) = 0
#
# Compile with:  See ex_05-setup.py for notes on installing and running
#
#                If you move this example to another directory, you must 
#                set the correct relative location for "braid.pyx" below.
#
# Help with:     $python3 ex_This is the simplest Cython example available, read the source
#
# Sample run:    From inside Python 3,     
#                >>> import ex_05
#                >>> core, app = ex_05.InitCoreApp()
#                >>> ex_05.run_Braid(core, app)
#
#                Output: 
#
#
# Sample 
# Parallel run:  See the file ex_05_run.py.  Run it from the shell with, 
#                $ mpirun -np 2  python3 ex_05_run.py
#
#                Output: 
#                Should be identical to serial run above
#                


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
    cdef object values

    def __cinit__(self, nx): 
        self.values = np.zeros((nx,), dtype='float') 


    cdef void SetVectorPtr(self, braid_Vector *u_ptr):
        ''' 
        This member function assigns the pointer value u_ptr to point to
        this vector. 
        '''
        u_ptr[0] = <braid_Vector> self


##
# Define user App to hold all "time-independent" information needed by the
# simulation.  This object is always available as the first parameter to any user
# function.
#
##
cdef class PyBraid_App:
    '''
    Use a Python class to define your braid_App
    '''

    cdef int rank      # MPI rank
    cdef int st        # time-stepper chosen (0: fwd Euler, 1: bwd Euler) 
    cdef int sc        # spatial coarsening used? (0: no,  1: yes) 
    cdef int nt        # number of points in space
    cdef double tstart # number of points in space
    cdef double tstop  # number of points in space
    cdef int nx        # number of points in space
    cdef double xL     # left initial condition in space and Dirichlet condition 
    cdef double xR     # right initial condition in space 
    cdef double eps    # diffusion coefficient
    cdef double a      # advection coefficient
    cdef object L      # spatial discretization matrix
    cdef object Phi    # list containing LinearOperators (for implicit) and sparse matrices 
                       # (for explicit) whose application allow for the evaluation of Phi 
    cdef object dts    # list containing dt value for each level 
    cdef object dxs    # list containing dx value for each level 
    
    def __cinit__(self, _rank): 
        self.rank = _rank


##
# Import Cython Braid Wrappers 
# (only do after declaration of Braid Vector and App)
include "../../braid/braid.pyx"


##
# Define user functions like Step, Init, ...
##

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    
    # Declare tstart and tstop as doubles, and then fill with values from Braid status
    tstart = 0.0 
    tstop = 0.0
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    dt = tstop - tstart

    # Declare level as an integer, then fill with value from Braid Status
    cdef int level = 0
    braid_StepStatusGetLevel(status, &level)

    # Cast u as a PyBraid_Vector, and do time-stepping.  Note
    # how the value array is written in-place
    pyU = <PyBraid_Vector> u

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app

    # Account for boundary values.  Note that matrix indexing pyApp.L[i,j] is
    # expensive.  Instead, we index directly into the the data array.  Also, 
    # only need right boundary condition if nonzero diffusion term.
    if pyApp.eps != 0.0:                                   
        pyU.values[0] += (dt*pyApp.L.data[2])*pyApp.xL     # Equivalent to   += (dt*pyApp.L[1,0])*pyApp.xL
        pyU.values[-1] += (dt*pyApp.L.data[-3])*pyApp.xR   # Equivalent to   += (dt*pyApp.L[-2,-1])*pyApp.xR
    else:
        pyU.values[0] += (dt*pyApp.L.data[1])*pyApp.xL     # Equivalent to   += (dt*pyApp.L[1,0])*pyApp.xL
    
    # Compute time-stepping matrix
    if pyApp.Phi[level] == []:
        I = eye(pyApp.nx, pyApp.nx, format='csr', dtype=pyApp.L.dtype)
        pyApp.Phi[level] = invert_sparse_mat_splu( I - dt*pyApp.L )
        pyApp.dts[level] = dt

    # Take step: make sure to write pyU in-place with [:]
    pyU.values[:] = pyApp.Phi[level] * pyU.values

    # Account for boundary conditions

    return 0


cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Allocate new vector
    pyU = PyBraid_Vector(pyApp.nx) 
    
    # Make sure to change values in-place with [:]
    if t == 0.0:
        nx = pyApp.nx
        pyU.values[0:(nx//10)] = pyApp.xL
        pyU.values[(nx//10):] = 0.0
        if pyApp.eps != 0.0:
            # Use right boundary condition if nonzero diffusion
            pyU.values[-1] = pyApp.xR

    else:
        pyU.values[:] = 6.123451
        # pyU.values[:] = np.random.rand(pyApp.nx)
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyU)
    
    # Set output, so that u_ptr points to pyU
    pyU.SetVectorPtr(u_ptr)

    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app

    # Allocate new vector
    pyV = PyBraid_Vector(pyApp.nx)
    
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Assign pyU's value to clone pyV
    #    Make sure to move values with the in-place operator [:]
    pyV.values[:] = pyU.values[:]
    
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
    
    # Compute AXPY, changing pyY in-place with [:]
    pyY.values[:] = alpha*pyX.values + beta*pyY.values

    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Compute norm 
    #  Note norm_ptr is a double array of size 1, and we index in at location [0]
    norm_ptr[0] = np.sqrt(np.dot(pyU.values, pyU.values))
    
    #print("Norm " + str(norm_ptr[0]) + ",  " + str(pyU.value[0]))
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app

    # Declare tindex as an integer, and the fill it in with the time index
    cdef int tindex
    braid_AccessStatusGetTIndex(status, &tindex)

    # Can get access to more information like
    # cdef double t
    # braid_AccessStatusGetT(status, &t)
    
    # Write this solution value to file in the format 
    #
    # ntime_steps
    # tstart
    # tstop
    # nspace_points
    # xstart
    # xstop
    # x[0]
    # x[1]
    #   .
    #   .
    #   .
    # x[k]
    filename = "%s.%05d.%04d"%("ex-05.out", tindex, pyApp.rank)
    f = open(filename, "w")
    f.write("%d\n"%pyApp.nt)
    f.write("%.14e\n"%pyApp.tstart)
    f.write("%.14e\n"%pyApp.tstop)
    f.write("%d\n"%pyApp.nx)
    f.write("%.14e\n"%0.0)
    f.write("%.14e\n"%1.0)

    for i in range(pyApp.nx):
        f.write( "%.14e\n"%pyU.values[i]);
    f.close()

    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    #  Note size_ptr is an integer array of size 1, and we index in at location [0]
    size_ptr[0] = pyApp.nx*sizeof(double)
    
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
    
    # Pack buffer, making sure to write in-place with [:]
    convert_carray_to_numpy(dbuffer, pyApp.nx)[:] = pyU.values[:]

    braid_BufferStatusSetSize(status, pyApp.nx*sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
   
    # Allocate new vector
    pyU = PyBraid_Vector(pyApp.nx)
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyU)
    
    # Unpack information in dbuffer, making sure to write in-place with [:]
    pyU.values[:] = convert_carray_to_numpy(dbuffer, pyApp.nx)[:]
    
    # Set output, so that u_ptr points to pyU
    pyU.SetVectorPtr(u_ptr)
    
    return 0

##
# Helper function to generate the spatial discretization for the
# advection-diffusion problem.  Different discretization types are supported.
def generate_spatial_disc(nx, xL, xR, a, eps, diff_discr='second_order', advect_discr='upwind'):
    
    # Assume unit interval with nx points
    dx = 1.0 / (nx - 1.0)

    # Construct diffusion matrix
    if diff_discr == 'second_order':
        # Construct the u_xx stencil of [1   -2   1]
        o = np.ones((3, nx))
        o[1,:] *= -2.0
        Diff =  (1./dx**2) * spdiags(o, np.array([1, 0,-1]), nx, nx).tocsr()
    else:
        exit("Advection discretization " + advect_discr + " unrecognized")
    
    # Construct advection matrix
    if advect_discr == 'upwind':
        # Construct the u_x stencil of [-1   1   0]
        o = np.ones((2, nx))
        o[0,:] *= -1.0
        
        Adv =  (1./dx) * spdiags(o, np.array([0,-1]), nx, nx).tocsr()
        
    elif advect_discr == 'central':
        # Construct the u_x stencil of [-1   0   1]
        o = np.ones((2, nx))
        o[0,:] *= -1.0
        
        Adv =  (1./(2.*dx)) * spdiags(o, np.array([1,-1]), nx, nx).tocsr()
    else:
        exit("Advection discretization " + advect_discr + " unrecognized")

    ##
    # Construct composite operator and return
    L = eps*Diff + a*Adv
    L.sort_indices()
    
    return L



# For a description of what each paramter here does, look at the help message
# printed by the driver, or look at the various braid_Set*() routines below. 
def InitCoreApp(print_help=False, ml=2, nu=1, nu0=1, CWt=1.0, skip=0, nx=16, ntime=60, xL=1.0,
        xR=0.0, eps=1.0, a=1.0, st=1, tol=1e-6, cf=2, mi=30, sc=0, fmg=0):

    cdef braid_Core core
    cdef double tstart
    cdef double tstop
    cdef MPI.Comm comm = MPI.COMM_WORLD
    cdef int rank
    
    if print_help:
        helpstring = "Braid simulation of 1D advection-diffusion\n\n" + \
             "Parameter options to InitCoreApp:\n" + \
             "  ml   <max_levels>   : set max levels\n" + \
             "  nu   <nrelax>       : set num F-C relaxations\n" + \
             "  nu0  <nrelax>       : set num F-C relaxations on level 0\n" + \
             "  CWt  <CWt>          : set C-relaxation weight on all levels\n" + \
             "  skip <set_skip>     : set skip relaxations on first down-cycle; 0: no skip;  1: skip\n" + \
             "  nx   <nspace>       : set num points in space\n" + \
             "  nt   <ntime>        : set num points in time\n" + \
             "  xL   <xLeft>        : set the left x-value (both as boundary condition and as initial condition)\n" + \
             "  xR   <xRight>       : set the right x-value (both as boundary condition and as initial condition)\n" + \
             "  eps  <epsilon>      : set the diffusion coefficient, u_t = a*u_x + eps*u_xx \n" + \
             "  a    <a>            : set the advection coefficient, u_t = a*u_x + eps*u_xx \n" + \
             "  st   <stepper>      : set the time stepper, 0: forward Euler, 1: backward Euler\n" + \
             "  tol  <tol>          : set stopping tolerance (scaled by sqrt(dt) sqrt(dx))\n" + \
             "  cf   <cfactor>      : set coarsening factor\n" + \
             "  mi   <max_iter>     : set max iterations\n" + \
             "  sc   <scoarsen>     : use spatial coarsening (bilinear) by factor of 2; must use 2^k sized grids\n" + \
             "  fmg                 : use FMG cycling \n"
        
        print(helpstring)

        return None,None
    
    else:
        # Set up braid core and app
        
        # number of time-steps, ntime, passed in as parameter
        tstart = 0.0
        tstop = 1.0
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
        
        # Set discretization options in App
        pyApp.st = st
        pyApp.sc = sc
        pyApp.nt = ntime
        pyApp.tstart = tstart
        pyApp.tstop = tstop
        pyApp.nx = nx
        pyApp.xL = xL
        pyApp.xR = xR
        pyApp.eps = eps
        pyApp.a = a
        pyApp.Phi = [ [] for i in range(ml) ]
        pyApp.dts = [ -1.0 for i in range(ml) ]
        
        # Define spatial discretization for each level
        pyApp.dxs = [] 
        for i in range(ml):
            pyApp.L = generate_spatial_disc(nx, xL, xR, a, eps)
            pyApp.dxs.append( (1.0 - 0.0) / (nx - 1.0) )

        # Scale tol by domain, assume x-domain is unit interval
        tol = tol/( np.sqrt((tstop - tstart)/(ntime-1.0))*np.sqrt((1.0 - 0.0)/(nx-1.0)) )
        
        # Set Braid options
        braid_SetMaxLevels(pyCore.getCore(), ml)
        braid_SetNRelax(pyCore.getCore(), -1, nu)
        braid_SetNRelax(pyCore.getCore(), 0, nu0)
        braid_SetCRelaxWt(pyCore.getCore(), -1, CWt)
        braid_SetSkip(pyCore.getCore(), skip)
        braid_SetAbsTol(pyCore.getCore(), tol)
        braid_SetCFactor(pyCore.getCore(), -1, cf)
        braid_SetAccessLevel(pyCore.getCore(), 1)
        braid_SetMaxIter(pyCore.getCore(), mi)
        if fmg == 1:
            braid_SetFMG(pyCore.getCore())
        
        return pyCore, pyApp


def run_Braid(PyBraid_Core pyCore, PyBraid_App pyApp):
    
    # Run Braid
    braid_Drive(pyCore.getCore())
    
    # Print per-level discretization information
    print("\n\n-----------------------------------------------------------------")
    print("-----------------------------------------------------------------\n")
    print("level       dx          dt        a*dt/dx      a*dt/dx^2"); 
    print("-----------------------------------------------------------------"); 
    for i in range(len(pyApp.Phi)): 
        dx = pyApp.dxs[i] 
        dt = pyApp.dts[i] 
        if dt == -1.0:
            break
        print(" %2d   |   %1.2e    %1.2e    %1.2e    %1.2e" %(i, dx, dt, pyApp.a*dt/dx, pyApp.a*dt/(dx*dx) )) 

    # Destroy App, decrementing smart pointer for app
    Py_DECREF(pyApp)
    del pyApp 

    # Destroy Braid
    braid_Destroy(pyCore.getCore())


