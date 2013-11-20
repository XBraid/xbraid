## Building Warp

-  The Warp examples assume that hypre and MFEM are installed in ../hypre and
   ../mfem, respectively.

-  To specify the compilers, flags and options for your machine, edit
   makefile.inc.  For now, we keep it simple and avoid using configure or
   cmake.

-  To make the library, libwarp.a,
   
         $ make

-  To make the examples
   
         $ make all

-  The makefile lets you pass some parameters like debug with 
   
         $ make debug=yes
   
   or
   
         $ make all debug=yes
   
   It would also be easy to add additional parameters, e.g., to compile with
   insure.  

-  To run the examples
   
         $ cd examples
         $ mpirun -np 4 drive-*

