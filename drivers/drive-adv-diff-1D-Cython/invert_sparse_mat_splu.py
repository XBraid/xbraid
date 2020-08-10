import scipy, numpy
from scipy import linalg
from scipy.sparse.linalg.interface import LinearOperator

__all__ = ['invert_sparse_mat_splu']


def invert_sparse_mat_splu(A, **kwargs):
    ''' 
    given a sparse matrix A, return a LinearOperator whose
    application is A-inverse by way of splu
    '''
    
    Acsc = A.tocsc()
    Acsc.eliminate_zeros()
    nonzero_cols = ( (Acsc.indptr[:-1] - Acsc.indptr[1:]) != 0).nonzero()[0]
    Map = scipy.sparse.eye(Acsc.shape[0], Acsc.shape[1], format='csc')
    Map = Map[:,nonzero_cols]
    Acsc = Map.T.tocsc()*Acsc*Map
    LU = scipy.sparse.linalg.splu(Acsc)
    LU_Map = Map

    def matvec(b):
        return LU_Map*LU.solve( numpy.ravel(LU_Map.T*b) )
    
    Ainv = LinearOperator(A.shape, matvec=matvec, dtype=A.dtype)

    return Ainv


