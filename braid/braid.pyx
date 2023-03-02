'''
Cython header file defining the Braid-Python interface
'''

cdef extern from "status.h":

    ##
    # Import the user-interface status structure definitions
    cdef struct _braid_Status_struct:
        pass
    ctypedef _braid_Status_struct *braid_Status

    cdef struct _braid_AccessStatus_struct:
        pass
    ctypedef _braid_AccessStatus_struct *braid_AccessStatus 

    cdef struct _braid_SyncStatus_struct:
        pass
    ctypedef _braid_SyncStatus_struct *braid_SyncStatus 

    cdef struct _braid_StepStatus_struct:
        pass
    ctypedef _braid_StepStatus_struct *braid_StepStatus

    cdef struct _braid_BufferStatus_struct:
        pass
    ctypedef _braid_BufferStatus_struct *braid_BufferStatus

    cdef struct _braid_CoarsenRefStatus_struct:
        pass
    ctypedef _braid_CoarsenRefStatus_struct *braid_CoarsenRefStatus
    

cdef extern from "braid_status.h":
    ##
    # Wrap StepStatus Routines
    int braid_StepStatusGetT (braid_StepStatus status, double *t_ptr)
    int braid_StepStatusGetTIndex (braid_StepStatus status, int *idx_ptr)
    int braid_StepStatusGetIter (braid_StepStatus status, int *iter_ptr)
    int braid_StepStatusGetLevel (braid_StepStatus status, int *level_ptr)
    int braid_StepStatusGetNLevels (braid_StepStatus status, int *nlevels_ptr)
    int braid_StepStatusGetNRefine (braid_StepStatus status, int *nrefine_ptr)
    int braid_StepStatusGetNTPoints (braid_StepStatus status, int *ntpoints_ptr)
    int braid_StepStatusGetTstop (braid_StepStatus status, double *tstop_ptr)
    int braid_StepStatusGetTstartTstop (braid_StepStatus status, double *tstart_ptr, double *tstop_ptr)
    int braid_StepStatusGetTol (braid_StepStatus status, double *tol_ptr)
    int braid_StepStatusGetRNorms (braid_StepStatus status, int *nrequest_ptr, double *rnorms_ptr)
    int braid_StepStatusGetOldFineTolx (braid_StepStatus status, double *old_fine_tolx_ptr)
    int braid_StepStatusGetSingleErrorEstStep (braid_StepStatus status, double *estimate)
    int braid_StepStatusGetCallingFunction (braid_StepStatus status, int *cfunction_ptr)
    int braid_StepStatusGetDeltaRank (braid_StepStatus status, int *rank_ptr)
    int braid_StepStatusGetBasisVec (braid_StepStatus status, braid_Vector *v_ptr, int index)
    int braid_StepStatusSetOldFineTolx (braid_StepStatus status, double old_fine_tolx)
    int braid_StepStatusSetTightFineTolx (braid_StepStatus status, double tight_fine_tolx)
    int braid_StepStatusSetRFactor (braid_StepStatus status, double rfactor)
    int braid_StepStatusSetRSpace (braid_StepStatus status, double r_space)

    ##
    # Wrap AccessStatus Routines
    int braid_AccessStatusGetT (braid_AccessStatus status, double *t_ptr) 
    int braid_AccessStatusGetTIndex (braid_AccessStatus status, int *idx_ptr)
    int braid_AccessStatusGetIter (braid_AccessStatus status, int *iter_ptr)
    int braid_AccessStatusGetLevel (braid_AccessStatus status, int *level_ptr)
    int braid_AccessStatusGetNLevels (braid_AccessStatus status, int *nlevels_ptr)
    int braid_AccessStatusGetNRefine (braid_AccessStatus status, int *nrefine_ptr)
    int braid_AccessStatusGetNTPoints (braid_AccessStatus status, int *ntpoints_ptr)
    int braid_AccessStatusGetResidual (braid_AccessStatus status, double *rnorm_ptr)
    int braid_AccessStatusGetDone (braid_AccessStatus status, int *done_ptr)
    int braid_AccessStatusGetTILD (braid_AccessStatus status, double *t_ptr, int *iter_ptr, int *level_ptr, int *done_ptr)
    int braid_AccessStatusGetWrapperTest (braid_AccessStatus status, int *wtest_ptr)
    int braid_AccessStatusGetCallingFunction (braid_AccessStatus status, int *cfunction_ptr)
    int braid_AccessStatusGetSingleErrorEstAccess (braid_AccessStatus status, double *estimate)
    int braid_AccessStatusGetDeltaRank (braid_AccessStatus status, int *rank_ptr)
    int braid_AccessStatusGetLocalLyapExponents (braid_AccessStatus status, double *exp_ptr, int *num_returned)
    int braid_AccessStatusGetBasisVec (braid_AccessStatus status, braid_Vector *v_ptr, int index)


    ##
    # Wrap CoarsenRefStatus Routines
    int braid_CoarsenRefStatusGetT (braid_CoarsenRefStatus status, double *t_ptr)
    int braid_CoarsenRefStatusGetTIndex (braid_CoarsenRefStatus status, int *idx_ptr)
    int braid_CoarsenRefStatusGetIter (braid_CoarsenRefStatus status, int *iter_ptr)
    int braid_CoarsenRefStatusGetLevel (braid_CoarsenRefStatus status, int *level_ptr)
    int braid_CoarsenRefStatusGetNLevels (braid_CoarsenRefStatus status, int *nlevels_ptr)
    int braid_CoarsenRefStatusGetNRefine (braid_CoarsenRefStatus status, int *nrefine_ptr)
    int braid_CoarsenRefStatusGetNTPoints (braid_CoarsenRefStatus status, int *ntpoints_ptr)
    int braid_CoarsenRefStatusGetCTprior (braid_CoarsenRefStatus status, double *ctprior_ptr)
    int braid_CoarsenRefStatusGetCTstop (braid_CoarsenRefStatus status, double *ctstop_ptr)
    int braid_CoarsenRefStatusGetFTprior (braid_CoarsenRefStatus status, double *ftprior_ptr)
    int braid_CoarsenRefStatusGetFTstop (braid_CoarsenRefStatus status, double *ftstop_ptr)
    int braid_CoarsenRefStatusGetTpriorTstop (braid_CoarsenRefStatus status, double *t_ptr, double *ftprior_ptr, double *ftstop_ptr, double *ctprior_ptr, double *ctstop_ptr)

    ##
    # Wrap BufferStatus Routines
    int braid_BufferStatusGetMessageType (braid_BufferStatus status, int *messagetype_ptr)
    int braid_BufferStatusGetTIndex ( braid_BufferStatus status, int *idx_ptr)
    int braid_BufferStatusGetLevel ( braid_BufferStatus status, int *level_ptr)
    int braid_BufferStatusSetSize ( braid_BufferStatus status, int size)
    int braid_BufferStatusSetBasisSize ( braid_BufferStatus status, int size)

    ##
    # Wrap SyncStatus Routines
    int braid_SyncStatusGetTIUL (braid_SyncStatus status, int *iloc_upper, int *iloc_lower, int level)
    int braid_SyncStatusGetTimeValues (braid_SyncStatus status, double **tvalues_ptr, int i_upper, int i_lower, int level)
    int braid_SyncStatusGetIter (braid_SyncStatus status, int *iter_ptr)
    int braid_SyncStatusGetLevel (braid_SyncStatus status, int *level_ptr)
    int braid_SyncStatusGetNLevels (braid_SyncStatus status, int *nlevels_ptr)
    int braid_SyncStatusGetNRefine (braid_SyncStatus status, int *nrefine_ptr)
    int braid_SyncStatusGetNTPoints (braid_SyncStatus status, int *ntpoints_ptr)
    int braid_SyncStatusGetDone (braid_SyncStatus status, int *done_ptr)
    int braid_SyncStatusGetCallingFunction (braid_SyncStatus status, int *cfunction_ptr)
    int braid_SyncStatusGetNumErrorEst (braid_SyncStatus status, int *npoints)
    int braid_SyncStatusGetAllErrorEst (braid_SyncStatus status, double *error_est)

cdef extern from "braid.h":
    
    ##
    # Wrap Braid Core
    cdef struct _braid_Core_struct:
        pass
    ctypedef _braid_Core_struct *braid_Core

    ## 
    # Wrap all the function pointers that User's can define
    ctypedef int (*braid_PtFcnStep)(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status)

    ctypedef int (*braid_PtFcnInit)(braid_App app, double t, braid_Vector *u_ptr)

    ctypedef int (*braid_PtFcnInitBasis)(braid_App app, double t, int index, braid_Vector *u_ptr)
    
    ctypedef int (*braid_PtFcnClone)(braid_App app, braid_Vector u, braid_Vector *v_ptr)

    ctypedef int (*braid_PtFcnFree)(braid_App app, braid_Vector u)

    ctypedef int (*braid_PtFcnSum)(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y)

    ctypedef int (*braid_PtFcnSpatialNorm)(braid_App app, braid_Vector u, double *norm_ptr)

    ctypedef int (*braid_PtFcnInnerProd)(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)

    ctypedef int (*braid_PtFcnAccess)(braid_App app, braid_Vector u, braid_AccessStatus status)

    ctypedef int (*braid_PtFcnBufSize)(braid_App app, int *size_ptr, braid_BufferStatus status)      

    ctypedef int (*braid_PtFcnBufPack)(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status)

    ctypedef int (*braid_PtFcnBufUnpack)(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status)

    ctypedef int (*braid_PtFcnResidual)(braid_App app, braid_Vector ustop, braid_Vector r, braid_StepStatus status) 

    ctypedef int (*braid_PtFcnSCoarsen)(braid_App app, braid_Vector fu, braid_Vector *cu_ptr, braid_CoarsenRefStatus  status)

    ctypedef int (*braid_PtFcnSRefine)(braid_App app, braid_Vector cu, braid_Vector *fu_ptr, braid_CoarsenRefStatus  status)

    ctypedef int (*braid_PtFcnSInit)(braid_App app, double t, braid_Vector *u_ptr)

    ctypedef int (*braid_PtFcnSClone)(braid_App app, braid_Vector u, braid_Vector *v_ptr)

    ctypedef int (*braid_PtFcnSFree)(braid_App app, braid_Vector  u)

    ctypedef int (*braid_PtFcnTimeGrid)(braid_App app, double *ta, int *ilower, int *iupper) 

    ##
    # Wrap BraidSet Routines
    int braid_PrintStats (braid_Core core)
    int braid_SetMaxLevels (braid_Core core, int max_levels)
    int braid_SetSkip (braid_Core core, int skip)
    int braid_SetRefine (braid_Core core, int refine)
    int braid_SetMaxRefinements (braid_Core core, int max_refinements)
    int braid_SetTPointsCutoff (braid_Core core, int tpoints_cutoff)
    int braid_SetMinCoarse (braid_Core core, int min_coarse)
    int braid_SetRelaxOnlyCG (braid_Core core, int relax_only_cg)
    int braid_SetAbsTol (braid_Core core, double atol)
    int braid_SetRelTol (braid_Core core, double rtol)
    int braid_SetNRelax (braid_Core core, int level, int nrelax)
    int braid_SetCRelaxWt(braid_Core  core, int level, double Cwt) 
    int braid_SetCFactor (braid_Core core, int level, int cfactor)
    int braid_SetMaxIter (braid_Core core, int max_iter)
    int braid_SetFMG (braid_Core core)
    int braid_SetNFMG (braid_Core core, int k)
    int braid_SetNFMGVcyc (braid_Core core, int nfmg_Vcyc)
    int braid_SetStorage (braid_Core core, int storage)
    int braid_SetTemporalNorm (braid_Core core, int tnorm)
    int braid_SetResidual (braid_Core core, braid_PtFcnResidual residual)
    int braid_SetFullRNormRes (braid_Core core, braid_PtFcnResidual residual)
    int braid_SetTimeGrid (braid_Core core, braid_PtFcnTimeGrid tgrid)
    int braid_SetSpatialCoarsen (braid_Core core, braid_PtFcnSCoarsen scoarsen)
    int braid_SetSpatialRefine (braid_Core core, braid_PtFcnSRefine srefine)
    int braid_SetPrintLevel (braid_Core core, int print_level)
    int braid_SetFileIOLevel (braid_Core core, int io_level)
    int braid_SetPrintFile (braid_Core core, const char *printfile_name)
    int braid_SetDefaultPrintFile (braid_Core core)
    int braid_SetAccessLevel (braid_Core core, int access_level)
    int braid_SplitCommworld (const libmpi.MPI_Comm *comm_world, int px, libmpi.MPI_Comm *comm_x, libmpi.MPI_Comm *comm_t)
    int braid_SetShell (braid_Core core, braid_PtFcnSInit sinit, braid_PtFcnSClone sclone, braid_PtFcnSFree sfree)
    int braid_SetRichardsonEstimation (braid_Core core, int est_error, int richardson, int local_order)
    int braid_SetDeltaCorrection (braid_Core core, int rank, braid_PtFcnInitBasis basis_init, braid_PtFcnInnerProd inner_prod)
    int braid_SetDeferDelta (braid_Core core, int level, int iter)
    int braid_SetLyapunovEstimation (braid_Core core, int relax, int cglv, int exponents)
    
    ##
    # Wrap BraidGet Routines
    int braid_GetNumIter (braid_Core core, int *iter_ptr)
    int braid_GetRNorms (braid_Core core, int *nrequest_ptr, double *rnorms)
    int braid_GetNLevels (braid_Core core, int *nlevels_ptr)
    int braid_GetSpatialAccuracy (braid_StepStatus status, double loose_tol, double tight_tol, double *tol_ptr)
    int braid_SetSeqSoln (braid_Core core, int seq_soln)
    int braid_GetMyID (braid_Core core, int *myid_ptr)

    ##
    # Wrap Braid Init, Drive, and Destroy
    int braid_Init (libmpi.MPI_Comm comm_world, libmpi.MPI_Comm comm, 
            double tstart, double tstop, int ntime, braid_App app, 
            braid_PtFcnStep step, braid_PtFcnInit init, braid_PtFcnClone clone, 
            braid_PtFcnFree free, braid_PtFcnSum sum, braid_PtFcnSpatialNorm spatialnorm, 
            braid_PtFcnAccess access,  braid_PtFcnBufSize bufsize, braid_PtFcnBufPack bufpack, 
            braid_PtFcnBufUnpack bufunpack, braid_Core *core_ptr)
    
    ##
    # Run Braid
    int braid_Drive (braid_Core core)
    
    ##
    # Deallocation
    int braid_Destroy (braid_Core core)


cdef object convert_carray_to_numpy(double * v, dim1, dim2=1, dim3=1):
    '''
    Helper function to cast C array v to an (dim1 x dim2 x dim3) 
    numpy array. Note, no copying of data is done.

    Input
    -----
    v     <double *>:   C data array to cast
    dim1  <int>     :   First dimension of v
    dim2  <int>     :   Second dimension of v (optional)
    dim3  <int>     :   Third dimension of v (optional)
    
    Output
    ----
    v_view          :   Numpy array that points to same memory as input v

    '''

    assert(dim1 > 0), "Array dimensions must be positive, dim1=%d"%dim1
    assert(dim2 >= 0), "Array dimensions must be positive, dim2=%d"%dim2
    assert(dim3 >= 0), "Array dimensions must be positive, dim3=%d"%dim3
    
    # Convert v to a linear array of size dim1*dim2*dim3
    dim = dim1*dim2*dim3
    cdef double[:] v_view = <double[:dim]> v
    v_np = np.asarray(v_view)

    if dim2 > 1 and dim3 == 1:
       v_np = v_np.reshape(dim1, dim2)        
    elif dim2 > 1 and dim3 > 1:
       v_np = v_np.reshape(dim1, dim2, dim3)        
    
    return v_np


cdef class PyBraid_Core:
    '''
    This class is a simple wrapper about the C-struct braid_Core
    This class allows us to return a braid_Core back into Python,
    and then have Python give that braid_Core back to Braid. 
    
    Note, that you cannot pass C-structs like braid_Core directly back to
    Python
    '''
    cdef braid_Core core

    def __init__(self): 
        pass

    cdef setCore(self, braid_Core _core):
        self.core = _core
    
    cdef braid_Core getCore(self):
        return self.core


cdef double* PyBraid_VoidToDoubleArray(void *buffer):
    '''
    Convert void * array to a double array 
    Note dbuffer is a C-array, so no bounds checking is done
    '''
    cdef double *dbuffer = <double*> buffer
    return dbuffer 


cdef int* PyBraid_VoidToIntArray(void *buffer):
    '''
    Convert void * array to a double array 
    Note ibuffer is a C-array, so no bounds checking is done
    '''
    cdef int *ibuffer = <int*> buffer
    return ibuffer 

