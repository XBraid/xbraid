#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#include "c_array.h"

/* definition of member functions */

/* create 1-d structures */
#undef CREATE
#define CREATE(type) \
type ## _array_1d *\
create_ ## type ## _array_1d( int n1 ){ \
  type ## _array_1d *tmp; \
  if (n1 <= 0){ \
    printf("ERROR in create_ ## type ## _array_1d: non-positive size!\n"); \
    return NULL; \
  } \
  if ((tmp = (type ## _array_1d *) malloc(sizeof(type ## _array_1d))) == NULL) \
    printf("memory error in type ## _array_1d\n"); \
  tmp->n1 = n1; \
  if ((tmp->arrayptr = (type *) malloc(n1*sizeof(type))) == NULL) \
    printf("memory error in type ## _array_1d, malloc\n"); \
  return tmp; \
}
/* define instances of the CREATE macro. This will define the functions */
/* create_int_array_1d, create_float_array_1d, create_double_array_1d and */
/* create_real_array_1d */
CREATE(int)
CREATE(float)
CREATE(double)
#ifndef NO_REAL
CREATE(real)
#endif
#undef CREATE

/* create 2-d structures */
#undef CREATE
#define CREATE(type) \
type ## _array_2d * \
create_ ## type ## _array_2d( int n1, int n2 ){ \
  type ## _array_2d *tmp; \
  if (n1*n2 <= 0){ \
    printf("ERROR in create_" #type "_array_2d: non-positive size!\n"); \
    return NULL; \
  } \
  if ((tmp = (type ## _array_2d *) malloc(sizeof(type ## _array_2d))) == NULL) \
    printf("memory error in create_" #type "_array_2d\n"); \
  tmp->n1 = n1; \
  tmp->n2 = n2; \
  if ((tmp->arrayptr = (type *) malloc(n1*n2*sizeof(type))) == NULL) \
    printf("memory error in create_" #type "_array_2d, malloc\n"); \
  return tmp; \
}
/* define instances of the CREATE macro. This will define the functions */
/* create_int_array_2d, create_float_array_2d, create_double_array_2d and */
/* create_real_array_2d */
CREATE(int)
CREATE(float)
CREATE(double)
#ifndef NO_REAL
CREATE(real)
#endif
#undef CREATE

/* create 3-d structures */
#undef CREATE
#define CREATE(type) \
type ## _array_3d * \
create_ ## type ## _array_3d( int n1, int n2, int n3){ \
  type ## _array_3d *tmp; \
  if (n1*n2*n3 <= 0){ \
    printf("ERROR in create_" #type "_array_3d: non-positive size!\n"); \
    return NULL; \
  } \
  if ((tmp = (type ## _array_3d *) malloc(sizeof(type ## _array_3d))) == NULL) \
    printf("memory error in " #type "_array_3d\n"); \
  tmp->n1 = n1; \
  tmp->n2 = n2; \
  tmp->n3 = n3; \
  if ((tmp->arrayptr = (type *) malloc(n1*n2*n3*sizeof(type))) == NULL) \
    printf("memory error in " #type "_array_3d, malloc\n"); \
  return tmp; \
}
/* define instances of the CREATE macro. This will define the functions */
/* create_int_array_3d, create_float_array_3d, create_double_array_3d and */
/* create_real_array_3d */
CREATE(int)
CREATE(float)
CREATE(double)
#ifndef NO_REAL
CREATE(real)
#endif
#undef CREATE


/* create 4-d structures */
#undef CREATE
#define CREATE(type) \
type ## _array_4d * \
create_ ## type ## _array_4d( int n1, int n2, int n3, int n4){ \
  type ## _array_4d *tmp; \
  if (n1*n2*n3*n4 <= 0){ \
    printf("ERROR in create_" #type "_array_4d: non-positive size!\n"); \
    return NULL; \
  } \
  if ((tmp = (type ## _array_4d *) malloc(sizeof(type ## _array_4d))) == NULL) \
    printf("memory error in create_" #type "_array_4d\n"); \
  tmp->n1 = n1; \
  tmp->n2 = n2; \
  tmp->n3 = n3; \
  tmp->n4 = n4; \
  if ((tmp->arrayptr = (type *) malloc(n1*n2*n3*n4*sizeof(type))) == NULL) \
    printf("memory error in create_" #type "_array_4d, malloc\n"); \
  return tmp; \
}
/* define instances of the CREATE macro. This will define the functions */
/* create_int_array_4d, create_float_array_4d, create_double_array_4d and */
/* create_real_array_4d */
CREATE(int)
CREATE(float)
CREATE(double)
#ifndef NO_REAL
CREATE(real)
#endif
#undef CREATE


/* delete 1-d array structures*/
#undef DELETE
#define DELETE(type) \
type ## _array_1d * \
delete_ ## type ## _array_1d( type ## _array_1d *struct_ptr){ \
  if (struct_ptr != NULL){ \
    free(struct_ptr->arrayptr); \
    free(struct_ptr); \
  } \
  return NULL; \
}
DELETE(int)
DELETE(float)
DELETE(double)
#ifndef NO_REAL
DELETE(real)
#endif
#undef DELETE

/* delete 2-d array structures*/
#undef DELETE
#define DELETE(type) \
type ## _array_2d * \
delete_ ## type ## _array_2d( type ## _array_2d *struct_ptr){ \
  if (struct_ptr != NULL){ \
    free(struct_ptr->arrayptr); \
    free(struct_ptr); \
  } \
  return NULL; \
}
DELETE(int)
DELETE(float)
DELETE(double)
#ifndef NO_REAL
DELETE(real)
#endif
#undef DELETE

/* delete 3-d array structures*/
#undef DELETE
#define DELETE(type) \
type ## _array_3d * \
delete_ ## type ## _array_3d( type ## _array_3d *struct_ptr){ \
  if (struct_ptr != NULL){ \
    free(struct_ptr->arrayptr); \
    free(struct_ptr); \
  } \
  return NULL; \
}
DELETE(int)
DELETE(float)
DELETE(double)
#ifndef NO_REAL
DELETE(real)
#endif
#undef DELETE

/* delete 4-d array structures*/
#undef DELETE
#define DELETE(type) \
type ## _array_4d * \
delete_ ## type ## _array_4d( type ## _array_4d *struct_ptr){ \
  if (struct_ptr != NULL){ \
    free(struct_ptr->arrayptr); \
    free(struct_ptr); \
  } \
  return NULL; \
}
DELETE(int)
DELETE(float)
DELETE(double)
#ifndef NO_REAL
DELETE(real)
#endif
#undef DELETE

/* output integer arrays to a binary file */

void write_int_array_1d( int_array_1d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, struct_ptr->n1 * sizeof(int));
}

void write_int_array_2d( int_array_2d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * sizeof(int));
}

void write_int_array_3d( int_array_3d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
  write( fd, (char *) &(struct_ptr->n3), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * struct_ptr->n3 * sizeof(int));
}

void write_int_array_4d( int_array_4d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
  write( fd, (char *) &(struct_ptr->n3), sizeof(int));
  write( fd, (char *) &(struct_ptr->n4), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * struct_ptr->n3 * struct_ptr->n4 * 
	 sizeof(int));
}

/* save integer arrays on ascii file */

void ascii_int_array_1d( int_array_1d *struct_ptr, FILE *fp, char *name){
  int i;
  fprintf( fp, "@int_array_1d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&int_array_1d_data\n" );
  for (i=1; i<= struct_ptr->n1; i++)
    fprintf( fp, "%i\n", compute_index_1d(struct_ptr, i) );
}

void ascii_int_array_2d( int_array_2d *struct_ptr, FILE *fp, char *name){
  int i, j;
  fprintf( fp, "@int_array_2d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&int_array_2d_data\n" );
  for (j=1; j <= struct_ptr->n2; j++)
    for (i=1; i <= struct_ptr->n1; i++)
      fprintf( fp, "%i\n", compute_index_2d(struct_ptr, i, j) );
}

void ascii_int_array_3d( int_array_3d *struct_ptr, FILE *fp, char *name){
  int i, j, k;
  fprintf( fp, "@int_array_3d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
  fprintf( fp, "&dim3 %i\n", struct_ptr->n3 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&int_array_3d_data\n" );
  for (k=1; k <= struct_ptr->n3; k++)
    for (j=1; j <= struct_ptr->n2; j++)
      for (i=1; i <= struct_ptr->n1; i++)
	fprintf( fp, "%i\n", compute_index_3d(struct_ptr, i, j, k) );
}

void ascii_int_array_4d( int_array_4d *struct_ptr, FILE *fp, char *name){
  int i, j, k, l;
  fprintf( fp, "@int_array_4d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
  fprintf( fp, "&dim3 %i\n", struct_ptr->n3 );
  fprintf( fp, "&dim4 %i\n", struct_ptr->n4 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&int_array_4d_data\n" );
  for (l=1; l <= struct_ptr->n4; l++)
    for (k=1; k <= struct_ptr->n3; k++)
      for (j=1; j <= struct_ptr->n2; j++)
	for (i=1; i <= struct_ptr->n1; i++)
	  fprintf( fp, "%i\n", compute_index_4d(struct_ptr, i, j, k, l) );
}

/* save real arrays on ascii file */

#ifndef NO_REAL
void ascii_real_array_1d( real_array_1d *struct_ptr, FILE *fp, char *name){
  int i;
  fprintf( fp, "@real_array_1d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&real_array_1d_data\n" );
  for (i=1; i<= struct_ptr->n1; i++)
    fprintf( fp, "%.18G\n", compute_index_1d(struct_ptr, i) );
}

void ascii_real_array_2d( real_array_2d *struct_ptr, FILE *fp, char *name){
  int i, j;
  fprintf( fp, "@real_array_2d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&real_array_2d_data\n" );
  for (j=1; j <= struct_ptr->n2; j++)
    for (i=1; i <= struct_ptr->n1; i++)
      fprintf( fp, "%.18G\n", compute_index_2d(struct_ptr, i, j) );
}

void ascii_real_array_3d( real_array_3d *struct_ptr, FILE *fp, char *name){
  int i, j, k;
  fprintf( fp, "@real_array_3d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
  fprintf( fp, "&dim3 %i\n", struct_ptr->n3 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&real_array_3d_data\n" );
  for (k=1; k <= struct_ptr->n3; k++)
    for (j=1; j <= struct_ptr->n2; j++)
      for (i=1; i <= struct_ptr->n1; i++)
	fprintf( fp, "%.18G\n", compute_index_3d(struct_ptr, i, j, k) );
}

void ascii_real_array_4d( real_array_4d *struct_ptr, FILE *fp, char *name){
  int i, j, k, l;
  fprintf( fp, "@real_array_4d %s\n", name );
/* output the size */
  fprintf( fp, "&dim1 %i\n", struct_ptr->n1 );
  fprintf( fp, "&dim2 %i\n", struct_ptr->n2 );
  fprintf( fp, "&dim3 %i\n", struct_ptr->n3 );
  fprintf( fp, "&dim4 %i\n", struct_ptr->n4 );
/* output the data */
  fprintf( fp, "# The first index changes the fastest\n" );
  fprintf( fp, "&real_array_4d_data\n" );
  for (l=1; l <= struct_ptr->n4; l++)
    for (k=1; k <= struct_ptr->n3; k++)
      for (j=1; j <= struct_ptr->n2; j++)
	for (i=1; i <= struct_ptr->n1; i++)
	  fprintf( fp, "%.18G\n", compute_index_4d(struct_ptr, i, j, k, l) );
}
#endif

/* dump array on stdout */
void 
print_int_array_1d( int_array_1d *struct_ptr ){
  int i;
/* output the size */
  printf( "Array dimensions: (%i)\n", struct_ptr->n1 );
/* output the data */
  printf( "Array contents:\n" );
  for (i=1; i <= struct_ptr->n1; i++)
    printf( "a(%i) = %i\n", i, compute_index_1d(struct_ptr, i) );
}

#ifndef NO_REAL
void 
print_real_array_1d( real_array_1d *struct_ptr ){
  int i;
/* output the size */
  printf( "Array dimensions: (%i)\n", struct_ptr->n1 );
/* output the data */
  printf( "Array contents:\n" );
  for (i=1; i <= struct_ptr->n1; i++)
    printf( "a(%i) = %.18G\n", i, compute_index_1d(struct_ptr, i) );
}
#endif

void 
print_int_array_2d( int_array_2d *struct_ptr ){
  int i, j;
/* output the size */
  printf( "Array dimensions: (%i,%i)\n", struct_ptr->n1, struct_ptr->n2 );
/* output the data */
  printf( "Array contents (first index is the fastest:\n" );
  for (j=1; j <= struct_ptr->n2; j++)
    {
      for (i=1; i <= struct_ptr->n1; i++)
	printf( "%i ", compute_index_2d(struct_ptr, i, j) );
      printf("\n");
    }
}

#undef PRINT_ARRAY
#define PRINT_ARRAY(type) \
void \
print_ ## type ## _array_2d( type ## _array_2d *struct_ptr ){ \
  int i, j; \
/* output the size */ \
  printf( "Array dimensions: (%i,%i)\n", struct_ptr->n1, struct_ptr->n2 ); \
/* output the data */ \
  printf( "Array contents:\n" ); \
  for (j=1; j <= struct_ptr->n2; j++) \
    for (i=1; i <= struct_ptr->n1; i++) \
      printf( "a(%i, %i) = %.18G\n", i, j, compute_index_2d(struct_ptr, i, j) ); \
}
#ifndef NO_REAL
PRINT_ARRAY(real)
#endif
PRINT_ARRAY(float)
PRINT_ARRAY(double)
#undef PRINT_ARRAY

void 
print_int_array_3d( int_array_3d *struct_ptr ){
  int i, j, k;
/* output the size */
  printf( "Array dimensions: (%i,%i,%i)\n", struct_ptr->n1, struct_ptr->n2, 
	 struct_ptr->n3 );
/* output the data */
  printf( "Array contents:\n" );
  for (k=1; k <= struct_ptr->n3; k++)
    for (j=1; j <= struct_ptr->n2; j++)
      for (i=1; i <= struct_ptr->n1; i++)
	printf( "a(%i, %i, %i) = %i\n", i, j, k, compute_index_3d(struct_ptr, i, j, k) );
}

#ifndef NO_REAL
void 
print_real_2d_element( real_array_2d *struct_ptr, int i, int j ){
/* output the data */
  printf( "Array (%i, %i) contents:", i, j );
  if (i >= 1 && i <= struct_ptr->n1 &&
      j >= 1 && j <= struct_ptr->n2)
    printf("%e\n", compute_index_2d(struct_ptr, i, j));
  else
    printf("out of bounds (1<=i<=%i, 1<=j<=%i)\n", struct_ptr->n1, struct_ptr->n2);
}

void 
print_real_3d_element( real_array_3d *struct_ptr, int i, int j, int k ){
/* output the data */
  printf( "Array (%i, %i, %i) contents:", i, j, k );
  if (i >= 1 && i <= struct_ptr->n1 &&
      j >= 1 && j <= struct_ptr->n2 &&
      k >= 1 && k <= struct_ptr->n3)
    printf("%e\n", compute_index_3d(struct_ptr, i, j, k));
  else
    printf("out of bounds (1<=i<=%i, 1<=j<=%i, 1<=k<=%i)\n", 
	   struct_ptr->n1, struct_ptr->n2, struct_ptr->n3);
}
#endif

void 
print_int_3d_element( int_array_3d *struct_ptr, int i, int j, int k ){
/* output the data */
  printf( "Array (%i, %i, %i) contents:", i, j, k );
  if (i >= 1 && i <= struct_ptr->n1 &&
      j >= 1 && j <= struct_ptr->n2 &&
      k >= 1 && k <= struct_ptr->n3)
    printf("%i\n", compute_index_3d(struct_ptr, i, j, k));
  else
    printf("out of bounds (1<=i<=%i, 1<=j<=%i, 1<=k<=%i)\n", 
	   struct_ptr->n1, struct_ptr->n2, struct_ptr->n3);
}


/* output real arrays on a binary file */

#ifndef NO_REAL
void write_real_array_1d( real_array_1d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, struct_ptr->n1 * sizeof(real));
}

void write_real_array_2d( real_array_2d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * sizeof(real));
}

void write_real_array_3d( real_array_3d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
  write( fd, (char *) &(struct_ptr->n3), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * struct_ptr->n3 * sizeof(real));
}

void write_real_array_4d( real_array_4d *struct_ptr, int fd){
/* output the size */
  write( fd, (char *) &(struct_ptr->n1), sizeof(int));
  write( fd, (char *) &(struct_ptr->n2), sizeof(int));
  write( fd, (char *) &(struct_ptr->n3), sizeof(int));
  write( fd, (char *) &(struct_ptr->n4), sizeof(int));
/* output the data */
  write( fd, (char *) struct_ptr->arrayptr, 
	 struct_ptr->n1 * struct_ptr->n2 * struct_ptr->n3 * struct_ptr->n4 * 
	 sizeof(real));
}
#endif

/* read an integer array from a binary file */

int_array_1d *read_int_array_1d( int fd ){
  int n1;
  int_array_1d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
/* create the array */
  tmp = create_int_array_1d( n1 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1*sizeof(int) );
/* return the pointer to the array */
  return tmp;
}

int_array_2d *read_int_array_2d( int fd ){
  int n1,n2;
  int_array_2d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
/* create the array */
  tmp = create_int_array_2d( n1, n2 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * sizeof(int) );
/* return the pointer to the array */
  return tmp;
}

int_array_3d *read_int_array_3d( int fd ){
  int n1, n2, n3;
  int_array_3d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
  read( fd, (char *) &n3, sizeof(int) );
/* create the array */
  tmp = create_int_array_3d( n1, n2, n3 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * n3 * sizeof(int) );
/* return the pointer to the array */
  return tmp;
}

int_array_4d *read_int_array_4d( int fd ){
  int n1, n2, n3, n4;
  int_array_4d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
  read( fd, (char *) &n3, sizeof(int) );
  read( fd, (char *) &n4, sizeof(int) );
/* create the array */
  tmp = create_int_array_4d( n1, n2, n3, n4 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * n3 * n4 * sizeof(int) );
/* return the pointer to the array */
  return tmp;
}


/* read a real array from a binary file */

#ifndef NO_REAL
real_array_1d *read_real_array_1d( int fd ){
  int n1;
  real_array_1d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
/* create the array */
  tmp = create_real_array_1d( n1 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1*sizeof(real) );
/* return the pointer to the array */
  return tmp;
}

real_array_2d *read_real_array_2d( int fd ){
  int n1,n2;
  real_array_2d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
/* create the array */
  tmp = create_real_array_2d( n1, n2 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * sizeof(real) );
/* return the pointer to the array */
  return tmp;
}

real_array_3d *read_real_array_3d( int fd ){
  int n1, n2, n3;
  real_array_3d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
  read( fd, (char *) &n3, sizeof(int) );
/* create the array */
  tmp = create_real_array_3d( n1, n2, n3 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * n3 * sizeof(real) );
/* return the pointer to the array */
  return tmp;
}

real_array_4d *read_real_array_4d( int fd ){
  int n1, n2, n3, n4;
  real_array_4d *tmp;
/* read the size */
  read( fd, (char *) &n1, sizeof(int) );
  read( fd, (char *) &n2, sizeof(int) );
  read( fd, (char *) &n3, sizeof(int) );
  read( fd, (char *) &n4, sizeof(int) );
/* create the array */
  tmp = create_real_array_4d( n1, n2, n3, n4 );
/* read the data */
  read( fd, (char *) tmp->arrayptr, n1 * n2 * n3 * n4 * sizeof(real) );
/* return the pointer to the array */
  return tmp;
}
#endif
