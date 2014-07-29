/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#ifndef C_ARRAY_H
#define C_ARRAY_H

#ifndef NO_REAL
#define NO_REAL
#endif

#include <stdio.h>

/* define real */
/*#include "real.h" */

/* note that all indexes are base 1, using fortran ordering for multi-dimensional arrays */
#define compute_index_1d(p,i) \
p->arrayptr[i-1]

#define compute_index_2d(p,i,j) \
p->arrayptr[i-1 + (j-1)*(p->n1)]

#define compute_index_3d(p,i,j,k) \
p->arrayptr[(i-1) + (j-1)*(p->n1) + (k-1)*(p->n1)*(p->n2)]

#define compute_index_4d(p,i,j,k,l) \
p->arrayptr[(i-1) + (j-1)*(p->n1) + (k-1)*(p->n1)*(p->n2) + \
	    (l-1)*(p->n1)*(p->n2)*(p->n3)]


/* 1-d structures */
#undef ARRAY
#define ARRAY(type) \
typedef struct { \
  int n1; \
  type *arrayptr; \
} type ## _array_1d;
/* now declare instances of the ARRAY macro. This will define int_array_1d, */
/* float_array_1d, double_array_1d and real_array_1d */
ARRAY(int)
ARRAY(float)
ARRAY(double)
#ifndef NO_REAL
ARRAY(real)
#endif
#undef ARRAY

/* 2-d structures */
#undef ARRAY
#define ARRAY(type) \
typedef struct { \
  int n1,n2; \
  type *arrayptr; \
} type ## _array_2d;
/* now declare instances of the ARRAY macro. This will define int_array_2d, */
/* float_array_2d, double_array_2d and real_array_2d */
ARRAY(int)
ARRAY(float)
ARRAY(double)
#ifndef NO_REAL
ARRAY(real)
#endif
#undef ARRAY

/* 3-d structures */
#undef ARRAY
#define ARRAY(type) \
typedef struct { \
  int n1,n2,n3; \
  type *arrayptr; \
} type ## _array_3d; 
/* now declare instances of the ARRAY macro. This will define int_array_3d, */
/* float_array_3d, double_array_3d and real_array_3d */
ARRAY(int)
ARRAY(float)
ARRAY(double)
#ifndef NO_REAL
ARRAY(real)
#endif
#undef ARRAY

/* 4-d structures */
#undef ARRAY
#define ARRAY(type) \
typedef struct { \
  int n1,n2,n3,n4; \
  type *arrayptr; \
} type ## _array_4d; 
/* now declare instances of the ARRAY macro. This will define int_array_4d, */
/* float_array_4d, double_array_4d and real_array_4d */
ARRAY(int)
ARRAY(float)
ARRAY(double)
#ifndef NO_REAL
ARRAY(real)
#endif
#undef ARRAY


/* member function specifications */
int_array_1d *
create_int_array_1d( int n1 );
int_array_2d *
create_int_array_2d( int n1, int n2 );
int_array_3d *
create_int_array_3d( int n1, int n2, int n3);
int_array_4d *
create_int_array_4d( int n1, int n2, int n3, int n4);

#ifndef NO_REAL
real_array_1d *
create_real_array_1d( int n1 );
real_array_2d *
create_real_array_2d( int n1, int n2 );
real_array_3d *
create_real_array_3d( int n1, int n2, int n3);
real_array_4d *
create_real_array_4d( int n1, int n2, int n3, int n4);
#endif

float_array_1d *
create_float_array_1d( int n1 );
float_array_2d *
create_float_array_2d( int n1, int n2 );
float_array_3d *
create_float_array_3d( int n1, int n2, int n3);
float_array_4d *
create_float_array_4d( int n1, int n2, int n3, int n4);

double_array_1d *
create_double_array_1d( int n1 );
double_array_2d *
create_double_array_2d( int n1, int n2 );
double_array_3d *
create_double_array_3d( int n1, int n2, int n3);
double_array_4d *
create_double_array_4d( int n1, int n2, int n3, int n4);

int_array_1d *
delete_int_array_1d( int_array_1d *struct_ptr);
int_array_2d *
delete_int_array_2d( int_array_2d *struct_ptr);
int_array_3d *
delete_int_array_3d( int_array_3d *struct_ptr);
int_array_4d *
delete_int_array_4d( int_array_4d *struct_ptr);

#ifndef NO_REAL
real_array_1d *
delete_real_array_1d( real_array_1d *struct_ptr);
real_array_2d *
delete_real_array_2d( real_array_2d *struct_ptr);
real_array_3d *
delete_real_array_3d( real_array_3d *struct_ptr);
real_array_4d *
delete_real_array_4d( real_array_4d *struct_ptr);
#endif

float_array_1d *
delete_float_array_1d( float_array_1d *struct_ptr);
float_array_2d *
delete_float_array_2d( float_array_2d *struct_ptr);
float_array_3d *
delete_float_array_3d( float_array_3d *struct_ptr);
float_array_4d *
delete_float_array_4d( float_array_4d *struct_ptr);

double_array_1d *
delete_double_array_1d( double_array_1d *struct_ptr);
double_array_2d *
delete_double_array_2d( double_array_2d *struct_ptr);
double_array_3d *
delete_double_array_3d( double_array_3d *struct_ptr);
double_array_4d *
delete_double_array_4d( double_array_4d *struct_ptr);

void 
write_int_array_1d( int_array_1d *struct_ptr, int fd);
void 
write_int_array_2d( int_array_2d *struct_ptr, int fd);
void 
write_int_array_3d( int_array_3d *struct_ptr, int fd);
void 
write_int_array_4d( int_array_4d *struct_ptr, int fd);

#ifndef NO_REAL
void 
write_real_array_1d( real_array_1d *struct_ptr, int fd);
void 
write_real_array_2d( real_array_2d *struct_ptr, int fd);
void 
write_real_array_3d( real_array_3d *struct_ptr, int fd);
void 
write_real_array_4d( real_array_4d *struct_ptr, int fd);
#endif

int_array_1d *
read_int_array_1d( int fd );
int_array_2d *
read_int_array_2d( int fd );
int_array_3d *
read_int_array_3d( int fd );
int_array_4d *
read_int_array_4d( int fd );

#ifndef NO_REAL
real_array_1d *
read_real_array_1d( int fd );
real_array_2d *
read_real_array_2d( int fd );
real_array_3d *
read_real_array_3d( int fd );
real_array_4d *
read_real_array_4d( int fd );
#endif

void 
ascii_int_array_1d( int_array_1d *struct_ptr, FILE *fp, char *name);
void 
ascii_int_array_2d( int_array_2d *struct_ptr, FILE *fp, char *name);
void 
ascii_int_array_3d( int_array_3d *struct_ptr, FILE *fp, char *name);
void 
ascii_int_array_4d( int_array_4d *struct_ptr, FILE *fp, char *name);

#ifndef NO_REAL
void 
ascii_real_array_1d( real_array_1d *struct_ptr, FILE *fp, char *name);
void 
ascii_real_array_2d( real_array_2d *struct_ptr, FILE *fp, char *name);
void 
ascii_real_array_3d( real_array_3d *struct_ptr, FILE *fp, char *name);
void 
ascii_real_array_4d( real_array_4d *struct_ptr, FILE *fp, char *name);
#endif

void 
print_int_array_1d( int_array_1d *struct_ptr );
void 
print_int_array_2d( int_array_2d *struct_ptr );
void 
print_int_array_3d( int_array_3d *struct_ptr );
#ifndef NO_REAL
void 
print_real_array_1d( real_array_1d *struct_ptr );
void 
print_real_array_2d( real_array_2d *struct_ptr );
#endif
void 
print_float_array_2d( float_array_2d *struct_ptr );
void 
print_double_array_2d( double_array_2d *struct_ptr );
void 
print_int_3d_element( int_array_3d *struct_ptr, int i, int j, int k );
#ifndef NO_REAL
void 
print_real_2d_element( real_array_2d *struct_ptr, int i, int j );
void 
print_real_3d_element( real_array_3d *struct_ptr, int i, int j, int k );
#endif

#endif
