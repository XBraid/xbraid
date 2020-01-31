/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision$
 ***********************************************************************EHEADER*/

/** \file mpistubs.h
 * \brief Fake MPI stubs to generate serial codes without MPI.
 *
 * This file contains the stubs to allow the user to fake using MPI, i.e., 
 * compile XBraid without MPI support for a purely serial version. 
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef braid_SEQUENTIAL

/******************************************************************************
 * MPI stubs to generate serial codes without mpi
 *****************************************************************************/

/* These types have associated creation and destruction routines */
//typedef int MPI_Comm;
typedef int MPI_Group;
typedef int MPI_Request;
typedef int MPI_Datatype;

typedef struct
{
   int MPI_SOURCE;
   int MPI_TAG;
} MPI_Status;
typedef int  MPI_Op;
typedef int  MPI_Aint;

#define  MPI_COMM_WORLD 0
#define  MPI_COMM_NULL  -1

#define  MPI_BOTTOM  0x0

#define  MPI_DOUBLE 0
#define  MPI_INT 1
#define  MPI_CHAR 2
#define  MPI_LONG 3
#define  MPI_BYTE 4
#define  MPI_REAL 5
#define  MPI_COMPLEX 6
#define  MPI_FLOAT 7

#define  MPI_SUM 0
#define  MPI_MIN 1
#define  MPI_MAX 2
#define  MPI_LOR 3
#define  MPI_SUCCESS 0
#define  MPI_STATUS_IGNORE 0

#define  MPI_UNDEFINED -9999
#define  MPI_REQUEST_NULL  0
#define  MPI_ANY_SOURCE    1
#define  MPI_ANY_TAG       1

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/* mpistubs.c */
int MPI_Init( int *argc , char ***argv );
int MPI_Finalize( void );
int MPI_Abort( MPI_Comm comm , int errorcode );
double MPI_Wtime( void );
double MPI_Wtick( void );
int MPI_Barrier( MPI_Comm comm );
int MPI_Comm_create( MPI_Comm comm , MPI_Group group , MPI_Comm *newcomm );
int MPI_Comm_dup( MPI_Comm comm , MPI_Comm *newcomm );
MPI_Comm MPI_Comm_f2c( int comm );
MPI_Comm MPI_Comm_c2f( int comm );
int MPI_Comm_size( MPI_Comm comm , int *size );
int MPI_Comm_rank( MPI_Comm comm , int *rank );
int MPI_Comm_free( MPI_Comm *comm );
int MPI_Comm_group( MPI_Comm comm , MPI_Group *group );
int MPI_Comm_split( MPI_Comm comm, int n, int m, MPI_Comm * comms );
int MPI_Group_incl( MPI_Group group , int n , int *ranks , MPI_Group *newgroup );
int MPI_Group_free( MPI_Group *group );
int MPI_Address( void *location , MPI_Aint *address );
int MPI_Get_count( MPI_Status *status , MPI_Datatype datatype , int *count );
int MPI_Alltoall( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int recvcount , MPI_Datatype recvtype , MPI_Comm comm );
int MPI_Allgather( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int recvcount , MPI_Datatype recvtype , MPI_Comm comm );
int MPI_Allgatherv( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int *recvcounts , int *displs , MPI_Datatype recvtype , MPI_Comm comm );
int MPI_Gather( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int recvcount , MPI_Datatype recvtype , int root , MPI_Comm comm );
int MPI_Gatherv( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int *recvcounts , int *displs , MPI_Datatype recvtype , int root , MPI_Comm comm );
int MPI_Scatter( void *sendbuf , int sendcount , MPI_Datatype sendtype , void *recvbuf , int recvcount , MPI_Datatype recvtype , int root , MPI_Comm comm );
int MPI_Scatterv( void *sendbuf , int *sendcounts , int *displs, MPI_Datatype sendtype , void *recvbuf , int recvcount , MPI_Datatype recvtype , int root , MPI_Comm comm );
int MPI_Bcast( void *buffer , int count , MPI_Datatype datatype , int root , MPI_Comm comm );
int MPI_Send( void *buf , int count , MPI_Datatype datatype , int dest , int tag , MPI_Comm comm );
int MPI_Recv( void *buf , int count , MPI_Datatype datatype , int source , int tag , MPI_Comm comm , MPI_Status *status );
int MPI_Isend( void *buf , int count , MPI_Datatype datatype , int dest , int tag , MPI_Comm comm , MPI_Request *request );
int MPI_Irecv( void *buf , int count , MPI_Datatype datatype , int source , int tag , MPI_Comm comm , MPI_Request *request );
int MPI_Send_init( void *buf , int count , MPI_Datatype datatype , int dest , int tag , MPI_Comm comm , MPI_Request *request );
int MPI_Recv_init( void *buf , int count , MPI_Datatype datatype , int dest , int tag , MPI_Comm comm , MPI_Request *request );
int MPI_Irsend( void *buf , int count , MPI_Datatype datatype , int dest , int tag , MPI_Comm comm , MPI_Request *request );
int MPI_Startall( int count , MPI_Request *array_of_requests );
int MPI_Probe( int source , int tag , MPI_Comm comm , MPI_Status *status );
int MPI_Iprobe( int source , int tag , MPI_Comm comm , int *flag , MPI_Status *status );
int MPI_Test( MPI_Request *request , int *flag , MPI_Status *status );
int MPI_Testall( int count , MPI_Request *array_of_requests , int *flag , MPI_Status *array_of_statuses );
int MPI_Wait( MPI_Request *request , MPI_Status *status );
int MPI_Waitall( int count , MPI_Request *array_of_requests , MPI_Status *array_of_statuses );
int MPI_Waitany( int count , MPI_Request *array_of_requests , int *index , MPI_Status *status );
int MPI_Allreduce( void *sendbuf , void *recvbuf , int count , MPI_Datatype datatype , MPI_Op op , MPI_Comm comm );
int MPI_Reduce( void *sendbuf , void *recvbuf , int count , MPI_Datatype datatype , MPI_Op op , int root , MPI_Comm comm );
int MPI_Scan( void *sendbuf , void *recvbuf , int count , MPI_Datatype datatype , MPI_Op op , MPI_Comm comm );
int MPI_Request_free( MPI_Request *request );
int MPI_Type_contiguous( int count , MPI_Datatype oldtype , MPI_Datatype *newtype );
int MPI_Type_vector( int count , int blocklength , int stride , MPI_Datatype oldtype , MPI_Datatype *newtype );
int MPI_Type_hvector( int count , int blocklength , MPI_Aint stride , MPI_Datatype oldtype , MPI_Datatype *newtype );
int MPI_Type_struct( int count , int *array_of_blocklengths , MPI_Aint *array_of_displacements , MPI_Datatype *array_of_types , MPI_Datatype *newtype );
int MPI_Type_commit( MPI_Datatype *datatype );
int MPI_Type_free( MPI_Datatype *datatype );

#endif
   
#ifdef __cplusplus
}
#endif

