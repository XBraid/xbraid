
#include <stdlib.h>
#include <stdio.h>

int
main( int argc, char *argv[] )
{
   int  arg_index;
   int  npoints, nprocs;
   int  proc, ilower, iupper, i, p, q, last_proc;

   if (argc > 2)
   {
      arg_index = 1;
      npoints = atoi(argv[arg_index++]);
      nprocs  = atoi(argv[arg_index++]);
   }
   else
   {
      printf("Problem with input\n");
      exit(1);
   }

   printf("\n npoints = %d, nprocs = %d\n", npoints, nprocs);

   /* GetDistribution */

   printf("\n GetDistribution:\n\n");
   for (proc = 0; proc < nprocs; proc++)
   {
      p = proc;
      ilower = npoints - ((nprocs-p)*npoints)/nprocs;
      p = proc+1;
      iupper = npoints - ((nprocs-p)*npoints)/nprocs - 1;

      printf(" %2d:", proc);
      for (i = ilower; i <= iupper; i++)
      {
         printf(" %2d", i);
      }
      printf("\n");
   }

   /* GetProc */

   printf("\n GetProc:\n");
   last_proc = -1;
   for (i = 0; i < npoints; i++)
   {
      proc = (i*nprocs)/npoints;
      
      if (proc > last_proc)
      {
         printf("\n");
         printf(" %2d:", proc);
      }
      printf(" %2d", i);

      last_proc = proc;
   }
   printf("\n\n");

   return (0);
}
