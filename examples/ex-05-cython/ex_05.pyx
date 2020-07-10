from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 
from scipy.sparse import spdiags, eye, linalg, csr_matrix, coo_matrix
from scipy.sparse.linalg.interface import LinearOperator
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
#                   u_t + a*u_x  =   eps*u_xx 
#                
#                on x in [0, 1], with t in [0, tstop=1.0]
#
#                Initial condition is a sine wave 
#
#                Boundary conditions are periodic 
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

    cdef int rank             # MPI rank
    cdef int sc               # spatial coarsening used? (0: no,  1: yes) 
    cdef int nt               # number of points in space
    cdef double tstart        # number of points in space
    cdef double tstop         # number of points in space
    cdef double eps           # diffusion coefficient
    cdef double a             # advection coefficient
    cdef object L             # list of spatial discretization matrices
    cdef object P             # list of spatial interpolation matrices 
    cdef object R             # list of spatial restriction matrices
    cdef object Phi           # list containing LinearOperators (for implicit) and sparse matrices 
                              # (for explicit) whose application allow for the evaluation of Phi 
    cdef object dts           # list containing dt value for each level 
    cdef object nx            # list containing the number of points in space for each level
    cdef object dxs           # list containing dx value for each level 
    cdef object time_discr    # string containing the desired time discretization
                              # 'BE', 'FE', 'SDIRK3', 'RK4'
    cdef object advect_discr  # string containing the desired advection discretization
                              # 'upwind', 'central', 'fourth', 'fourth_diss', 'fourth_diss_sq'
    cdef object diff_discr    # string containing the desired diffusion discretization
                              # 'second_order', 'fourth_order'
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

    L =  pyApp.L[level]
    nx = pyApp.nx[level]

    # Compute time-stepping matrix
    if pyApp.Phi[level] == []:
        pyApp.dts[level] = dt
        
        if pyApp.time_discr == 'BE':
            # Compute time-stepping operator for Backward Euler
            I = eye(nx, nx, format='csr', dtype=L.dtype)
            Inv = invert_sparse_mat_splu( I - dt*L )

            def matvec(y):
                return Inv*y #  if boundary conditions, Inv*(y + dt*pyApp.bcs)
        
            pyApp.Phi[level] = LinearOperator(L.shape, matvec=matvec, dtype=L.dtype)


        elif pyApp.time_discr == 'SDIRK3':
            # Compute time-stepping operator for SDIRK3 
            # 
            #   a  |   a
            #   c  |  c-a    a
            #   1  |   b   1-a-b  a
            # -----+----------------
            #      |   b   1-a-b  a
            a = 0.435866521508458999416019
            b = 1.20849664917601007033648
            c = 0.717933260754229499708010
            a_21 = c-a; a_31 = b; a_32 = 1-a-b
            b_1 = b; b_2 = 1-a-b; b_3 = a

            # Use SPLU to get our inverse 
            I = eye(nx, nx, format='csr', dtype=L.dtype)
            SDIRK_Diag_Inv = invert_sparse_mat_splu( I - dt*a*L )
            
            # Define a mat-vec that applies SDIRK3
            def matvec(y):
                
                k0 = y
                
                k1 = L*k0                             # if boundary conditions,  + pyApp.bcs
                k1 = SDIRK_Diag_Inv*( k1 )                  
                                                            
                k2 = L*(k0 + dt*a_21*k1)              # if boundary conditions,  + pyApp.bcs
                k2 = SDIRK_Diag_Inv*( k2 )
                
                k3 = L*(k0 + dt*a_31*k1 + dt*a_32*k2) # if boundary conditions,  + pyApp.bcs
                k3 = SDIRK_Diag_Inv*( k3 ) 
                
                return k0 + dt*(b_1*k1 + b_2*k2 + b_3*k3)
            
            pyApp.Phi[level] = LinearOperator(L.shape, matvec=matvec, dtype=L.dtype)

        elif pyApp.time_discr == 'FE':
            # Compute time-stepping operator for forward Euler
            I = eye(nx, nx, format='csr', dtype=L.dtype)
            Phi_Mat = I + dt*L 

            def matvec(y):
                return Phi_Mat*y #  if boundary conditions, Phi_Mat*(y + dt*pyApp.bcs)
        
            pyApp.Phi[level] = LinearOperator(L.shape, matvec=matvec, dtype=L.dtype)

        elif pyApp.time_discr == 'RK4':
            
            # Compute time-stepping operator for RK-4
            def matvec(y):
                k1 = L*y
                k2 = L*(y + dt/2*k1)
                k3 = L*(y + dt/2*k2)
                k4 = L*(y + dt*k3)
                y =  y + (dt/6.)*(k1 + 2*k2 + 2*k3 + k4)
            
                return y                    
            
            pyApp.Phi[level] = LinearOperator(L.shape, matvec=matvec, dtype=L.dtype)

    # Take step: make sure to write pyU in-place with [:]
    pyU.values[:] = pyApp.Phi[level] * pyU.values

    return 0


cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Allocate new vector
    nx = pyApp.nx[0]
    pyU = PyBraid_Vector(nx) 
    
    # Make sure to change values in-place with [:]
    if t == 0.0:
        mesh_x = np.linspace(0,1.0,nx)
        indys = np.array( (mesh_x >= 0.25), dtype=int) + np.array( (mesh_x <= 0.5), dtype=int) == 2
        pyU.values = np.sin(2*np.pi*mesh_x)
        #pyU.values[indys] = -np.sin(4*np.pi*mesh_x[indys])
        
        #pyU.values[:] = 0.0
        #pyU.values[(nx//10):(nx//5)] = 1.0
        ## If boundary conditions, need to set [0] and [-1] entries

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

    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    
    # Allocate new vector
    pyV = PyBraid_Vector(pyU.values.shape[0])
    
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
    f.write("%d\n"%pyApp.nx[0])
    f.write("%.14e\n"%0.0)
    f.write("%.14e\n"%1.0)

    for i in range(pyApp.nx[0]):
        f.write( "%.14e\n"%pyU.values[i]);
    f.close()

    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    #  Note size_ptr is an integer array of size 1, and we index in at location [0]
    size_ptr[0] = (pyApp.nx[0] + 1)*sizeof(double)
    
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    # Cast u as a PyBraid_Vector
    pyU = <PyBraid_Vector> u
    nx = pyU.values.shape[0]

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
    
    # Pack buffer, making sure to write in-place with [:]
    convert_carray_to_numpy(dbuffer, nx+1)[0] = nx 
    convert_carray_to_numpy(dbuffer, nx+1)[1:] = pyU.values[:]

    braid_BufferStatusSetSize(status, (nx+1)*sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    
    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # Convert void * to a double array (note dbuffer is a C-array, so no bounds checking is done) 
    dbuffer = PyBraid_VoidToDoubleArray(buffer)
   
    # Allocate new vector
    nx = int(convert_carray_to_numpy(dbuffer, 1)[0])
    pyU = PyBraid_Vector(nx)
    
    # Must increment smart pointer for pyU, or else Python will delete it 
    Py_INCREF(pyU)
    
    # Unpack information in dbuffer, making sure to write in-place with [:]
    pyU.values[:] = convert_carray_to_numpy(dbuffer, nx+1)[1:]
    
    # Set output, so that u_ptr points to pyU
    pyU.SetVectorPtr(u_ptr)
    
    return 0


## 1D Bilinear interpolation
## Input argument is nc, the number of points on the coarse grid
## nc must be a power of 2 plus 1
def interpolation1d(nc):
    d = np.repeat([[1, 2, 1]], nc, axis=0).T
    I = np.zeros((3,nc), dtype=int)
    for i in range(nc):
        I[:,i] = [2*i, 2*i+1, 2*i+2]
    J = np.repeat([np.arange(nc)], 3, axis=0)
    P = coo_matrix( (d.ravel(), (I.ravel(), J.ravel()))).tocsr()
    return 0.5 * P

## 1D Injection interpolation
## Input argument is nc, the number of points on the coarse grid
def injection1d(nc):
    d = np.repeat([[1]], nc, axis=0).T
    I = np.zeros((1,nc), dtype=int)
    for i in range(nc):
        I[:,i] = [2*i+1]
    J = np.repeat([np.arange(nc)], 1, axis=0)
    P = coo_matrix( (d.ravel(), (I.ravel(), J.ravel())), shape=(2*nc+1, nc)).tocsr()
    return P


cdef int my_coarsen(braid_App app, braid_Vector fu, braid_Vector *cu_ptr, braid_CoarsenRefStatus status):
    # Cast fu as a PyBraid_Vector
    pyFU = <PyBraid_Vector> fu

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # This returns the level for fu
    cdef int level = 0
    braid_CoarsenRefStatusGetLevel(status, &level)

    # Allocate new vector
    pyCU = PyBraid_Vector(pyApp.nx[level+1])
    
    # Restrict, making sure to write in-place
    pyCU.values[:] = pyApp.R[level]*pyFU.values
    
    # Must increment smart pointer for pyCU, or else Python will delete it 
    Py_INCREF(pyCU)
    
    # Set output, so that cu_ptr points to pyCU
    pyCU.SetVectorPtr(cu_ptr)


cdef int my_refine(braid_App app, braid_Vector cu, braid_Vector *fu_ptr, braid_CoarsenRefStatus status):
    # Cast cu as a PyBraid_Vector
    pyCU = <PyBraid_Vector> cu

    # Cast app as a PyBraid_App
    pyApp = <PyBraid_App> app
    
    # This returns the level for fu_ptr
    cdef int level = 0
    braid_CoarsenRefStatusGetLevel(status, &level)

    # Allocate new vector
    pyFU = PyBraid_Vector(pyApp.nx[level])
    
    # Interpolate, making sure to write in-place
    pyFU.values[:] = pyApp.P[level]*pyCU.values
    
    # Must increment smart pointer for pyFU, or else Python will delete it 
    Py_INCREF(pyFU)
    
    # Set output, so that fu_ptr points to pyFU
    pyFU.SetVectorPtr(fu_ptr)




##
# Helper function to generate the spatial discretization for the
# advection-diffusion problem.  Different discretization types are supported.
def generate_spatial_disc(nx, a, eps, diff_discr='second_order', advect_discr='upwind'):
    
    # Assume unit interval with nx points
    dx = 1.0 / (nx - 1.0)

    # Declare basic difference matrices.  # To avoid any boundary effects, we
    # add a padding of 10 (to be removed later)
    o = np.ones((2,nx+10))
    o[0,:] *= -1.
    #
    Dplus = (1./dx)*spdiags(o, [0,1], nx+10, nx+10, format='csr')
    Dminus =(1./dx)*spdiags(o, [-1,0], nx+10, nx+10, format='csr')
    Dzero = (1./(2.*dx))*spdiags(o, [-1,1], nx+10, nx+10, format='csr')  # equiv to Dzero = 0.5*(Dplus + Dminus)
    I = eye(nx+10, nx+10, format='csr')

    # Construct diffusion matrix
    if diff_discr == 'second_order':
        # Construct the u_xx stencil of [1   -2   1]
        Diff =  Dplus*Dminus
    elif diff_discr == 'fourth_order':
        # Construct the fourth order u_xx stencil 
        Diff =  Dplus*Dminus*(I - ((dx**2)/12.)*Dplus*Dminus)
    else:
        exit("Diffusion discretization " + diff_discr + " unrecognized")


    # Construct advection matrix
    if advect_discr == 'upwind':
        # Construct the u_x stencil of [-1   1   0]
        Adv =  Dminus 

    elif advect_discr == 'central':
        # Construct the u_x stencil of [-1   0   1]
        Adv =  Dzero

    elif (advect_discr == 'fourth') or (advect_discr == 'fourth_diss') or (advect_discr == 'fourth_diss_sq'):
       
        # Set fourth order centered stecil for d/dx
        Adv = Dzero*(I - ((dx**2)/6.)*Dplus*Dminus)
        
        # Add fourth-order hyperviscosity
        if advect_discr == 'fourth_diss':
            # Add one order of artificial diss (third order in space)
            Adv += 1.1*(dx**3)*Dplus*Dplus*Dminus*Dminus

        elif advect_discr == 'fourth_diss_sq':
            # Add two orders of artificial diss (second order in space)
            Adv += 1.1*(dx**2)*Dplus*Dplus*Dminus*Dminus
               
        # Note: if no periodicity in space, you need to fix up the last two
        # rows of Adv, so that they represent accurate first derivative
        # computations (because there is no BC to account for the truncated
        # stencils).  You could use one-sided stenicls, or just make the last
        # two rows upwinding.

    else:
        exit("Advection discretization " + advect_discr + " unrecognized")

    ##
    # Construct composite operator, remember that we solve
    #    u_t + a*u_x  =   eps*u_xx 
    # but we are forming here
    #    u_t = L u = -a u_x + eps u_xx
    L = - a*Adv + eps*Diff
    
    ##
    # Remove padding around boundary
    # Note, because we use periodic boundary conditions, the computations
    # of Dplus*Dminus (and so on) are inaccurate near the boundary.  Thus,
    # we inserted a plus 10 padding above, that we now remove.
    L = L[5:(5+nx), 5:(5+nx)]

    ##
    # Enforce periodic BCs
    L = L.tolil()
    L[0,-1] = L[2,1]
    L[0,-2] = L[2,0]
    L[1,-1] = L[2,0]
    #
    L[-1,0] = L[2,3]
    L[-1,1] = L[2,4]
    L[-2,0] = L[2,4]
    
    L = L.tocsr()
    L.sort_indices()
    return L


# For a description of what each paramter here does, look at the help message
# printed by the driver, or look at the various braid_Set*() routines below. 
def InitCoreApp(print_help=False, ml=2, nu=1, nu0=1, CWt=1.0, skip=0, nx=16, ntime=60, 
        eps=1.0, a=1.0, tol=1e-6, cf=2, mi=30, sc=0, fmg=0, advect_discr='upwind', 
        diff_discr='second_order', time_discr='BE'):

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
             "  eps  <epsilon>      : set the diffusion coefficient, u_t = a*u_x + eps*u_xx \n" + \
             "  a    <a>            : set the advection coefficient, u_t = a*u_x + eps*u_xx \n" + \
             "  tol  <tol>          : set stopping tolerance (scaled by sqrt(dt) sqrt(dx))\n" + \
             "  cf   <cfactor>      : set coarsening factor\n" + \
             "  mi   <max_iter>     : set max iterations\n" + \
             "  sc   <scoarsen>     : use spatial coarsening (bilinear) by factor of 2; must use 2^k sized grids\n" + \
             "  fmg                 : use FMG cycling \n" + \
             "  time_discr          : time discretization to use, choose one of \n" + \
             "                      :   'BE', 'FE', 'SDIRK3', 'RK4' \n" + \
             "  diff_discr          : diffusion discretization to use, choose one of \n" + \
             "                      :   'second_order', 'fourth_order' \n" + \
             "  advect_discr        : spatial discretization for advection term\n" +\
             "                      :   The advection term can be one of \n" + \
             "                      :   'upwind', 'central', 'fourth', 'fourth_diss', 'fourth_diss_sq' \n" + \
             "                      :   where these options are upwinding, central differencing, fourth order differencing, \n" + \
             "                      :   fourth order differencing with one order of artificial dissipation, and \n" + \
             "                      :   fourth order differencing with two orders of artificial dissipation.  \n\n" + \
             "                      :   The diffusion term is always discretized with classic second order differencing. \n"
        
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
        pyApp.sc = sc
        pyApp.nt = ntime
        pyApp.tstart = tstart
        pyApp.tstop = tstop
        pyApp.eps = eps
        pyApp.a = a
        pyApp.Phi = [ [] for i in range(ml) ]
        pyApp.dts = [ -1.0 for i in range(ml) ]
        pyApp.time_discr = time_discr
        pyApp.advect_discr = advect_discr
        pyApp.diff_discr = diff_discr
        
        # Define spatial discretization for each level
        pyApp.nx = [] 
        pyApp.dxs = [] 
        pyApp.L = []
        pyApp.R = []
        pyApp.P = []
        pyApp.nx.append(nx)
        for i in range(ml):
            pyApp.L.append(generate_spatial_disc(pyApp.nx[i], a, eps, advect_discr=pyApp.advect_discr, diff_discr=pyApp.diff_discr))
            
            if(pyApp.sc == 1):
                nc = (nx // 2**(i+1)) + 1
                pyApp.R.append(injection1d(nc))
                pyApp.P.append(interpolation1d(nc))
            else:
                nc = nx

            pyApp.dxs.append( (1.0 - 0.0) / (pyApp.nx[i] - 1.0) )
            if i != ml-1:
                pyApp.nx.append( nc )

        #  NNNeed some way to not coarsen beyond a grid of 3 points ... or just
        #leave that to the user?  or just stop after 3 points and do the
        #identity?  
        # Set spatial coarsening
        if(pyApp.sc == 1):
            # Check that the spatial grid size is a power of 2 plus 1
            if (2**np.floor(np.log2(pyApp.nx)) + 1) != pyApp.nx:
                exit("Spatial coarsening only works with power of 2 plus 1 grids")

            #braid_SetSpatialCoarsen(pyCore.getCore(), my_coarsen)
            #braid_SetSpatialRefine(pyCore.getCore(), my_refine)


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


