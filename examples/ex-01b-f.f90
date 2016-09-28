!BHEADER
!* Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
!* Produced at the Lawrence Livermore National Laboratory. Written by 
!* Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
!* Dobrev, et al. LLNL-CODE-660355. All rights reserved.
!* 
!* This file is part of XBraid. Email xbraid-support@llnl.gov for support.
!* 
!* This program is free software; you can redistribute it and/or modify it under
!* the terms of the GNU General Public License (as published by the Free Software
!* Foundation) version 2.1 dated February 1999.
!* 
!* This program is distributed in the hope that it will be useful, but WITHOUT ANY
!* WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
!* License for more details.
!* 
!* You should have received a copy of the GNU Lesser General Public License along
!* with this program; if not, write to the Free Software Foundation, Inc., 59
!* Temple Place, Suite 330, Boston, MA 02111-1307 USA
!*
!EHEADER


!  Example 01b (in Fortran 90)
!
!  Compile with: make ex-01b-f
!
!  Sample run:   mpirun -np 2 ex-01b-f
!
!  Description:
!
!  Solve the scalar ODE: u' = lambda u
!  
!  Same as ex-01, only show how to implement more advanced XBraid features.
!  
!  When run with the default 10 time steps, the solution is:
!  $ ./ex-01b-f
!  $ cat ex-01b-f.out.00*
!    1.00000000000000e+00
!    5.00000000000000e-01
!    2.50000000000000e-01
!    1.25000000000000e-01
!    6.25000000000000e-02
!    3.12500000000000e-02
!    1.56250000000000e-02
!    7.81250000000000e-03
!    3.90625000000000e-03
!    1.95312500000000e-03
!    9.76562500000000e-04


! F90 modules are a convenient way of defining XBraid vectors and app structure
module braid_types
   
   ! Declare vector object (just a scalar)
   type my_vector
      double precision val
   end type my_vector
   
   ! Declare app object
   ! You can put anything in app
   type my_app
      double precision              :: tstart
      double precision              :: tstop
      double precision, allocatable :: dt(:)
      integer                       :: mydt
      integer                       :: ntime
      integer                       :: comm
      integer                       :: bufsize
      integer                       :: pt
      integer                       :: rank
   end type my_app
   
   ! Many Fortran compilers have a sizeof( ) function but its not part of the
   ! language.  So, define variable sizes here.
   integer, parameter :: sizeof_double = 8
   integer, parameter :: sizeof_int = 4

   ! braid_Core will store a pointer value to the braid_Core structure
   ! makre sure that the correct kind() is chosen so that it will hold a memory
   ! address.  This holds true for all the other similar pointers below, such as
   ! the status structures
   integer (kind=8)  :: braid_core

end module braid_types


! Helper function: Replace character a with character b
subroutine replace(s, a, b, length)
   implicit none

   integer :: length, i
   character(len=length) s
   character(1) a, b

   do i = 1, length
      if( s(i:i) == a) then
         s(i:i) = b
      endif
   end do
end subroutine replace


subroutine init_timesteps(app)

   use braid_types
   implicit none
   type(my_app)     :: app
   integer          :: ntime, i
   double precision :: dt
   
   ntime = app%ntime
   dt    = (app%tstop - app%tstart) / app%ntime
   
   allocate(app%dt(ntime))

   ! example with varying time step size
   if (app%mydt == 2) then
      do i = 1,int(0.5*ntime)
         app%dt(i) = dt * dble(0.5);
      enddo
      do i = int(0.5*ntime)+1,ntime
         app%dt(i) = dt * dble(1.5);
      enddo
   ! default to constant time step size
   else
      do i = 1,ntime
         app%dt(i) = dt;
      enddo
   endif

end subroutine

subroutine print_timegrid(app, ta, ilower, iupper)

   use braid_types
   implicit none
   type(my_app)                             :: app
   double precision, dimension(app%ntime+1) :: ta
   integer                                  :: ilower, iupper
   integer                                  :: i, out_unit
   character(len=5)  :: rank_string
   character(len=25) :: fname = "timegrid"
   character(len=8)  :: fname_short = "timegrid"
   character(len=1)  :: dot = "."

   write(rank_string, "(I5)")  app%rank
   call replace(rank_string, ' ', '0', 5)
   fname = fname_short // dot // rank_string
   out_unit = 11
   open (unit=out_unit, file=fname, action="write", status="replace")
   do i = ilower,iupper
      write(out_unit,*) ta(i-ilower+1)
   enddo
   close (out_unit)

end subroutine


subroutine braid_timegrid_f90(app, ta, ilower, iupper)

   use braid_types
   implicit none
   type(my_app)                             :: app
   double precision, dimension(app%ntime+1) :: ta
   integer                                  :: ilower, iupper
   integer                                  :: i
   double precision                         :: tstart ! time corresponding to ilower, i.e. lower time index value for this processor

   ! Start from the global tstart to compute the local tstart
   tstart = app%tstart
   do i = 1,ilower
      tstart = tstart + app%dt(i)
   enddo
   ! Assign time point values for local time point index values ilower:iupper
   do i = ilower+1,iupper+1
      ta(i-ilower) = tstart
      tstart       = tstart + app%dt(i)
   enddo
   call print_timegrid(app, ta, ilower, iupper)

end subroutine


subroutine braid_Init_Vec_F90(app, t, u_ptr)
   
   ! Braid types
   use braid_types
   implicit none
   type(my_vector), pointer :: u_ptr
   type(my_app)             :: app
   
   ! Other declarations
   double precision :: t, val
   
   ! Initialize u_ptr
   allocate(u_ptr)
   if (t == 0.0) then
      u_ptr%val = 1.0
   else
      !call random_number(val)
      val = 0.456
      u_ptr%val = val
   endif

end subroutine braid_Init_Vec_F90


subroutine braid_Free_F90(app, u)
   
   ! Braid types
   use braid_types
   implicit none
   type(my_app)             :: app
   type(my_vector), pointer :: u
   
   deallocate(u) 
end subroutine braid_Free_F90


subroutine braid_Clone_F90(app, u, v_ptr)
   
   ! Braid types
   use braid_types
   implicit none
   type(my_vector)          :: u
   type(my_vector), pointer :: v_ptr
   type(my_app)             :: app
   
   ! Clone u into v_ptr
   allocate(v_ptr)
   v_ptr%val = u%val

end subroutine braid_Clone_F90


subroutine braid_Sum_F90(app, alpha, x, beta, y)
   
   ! Braid types
   use braid_types
   implicit none
   type(my_vector)          :: x, y
   type(my_app)             :: app
   
   ! Other declarations
   double precision alpha, beta
   
   y%val = alpha*(x%val) + beta*(y%val)

end subroutine braid_Sum_F90


subroutine braid_SpatialNorm_F90(app, u, norm_ptr)
   
   ! Braid types
   use braid_types
   implicit none
   type(my_vector)          :: u
   type(my_app)             :: app
   
   ! Other declarations
   double precision norm_ptr
 
   norm_ptr = sqrt( (u%val)*(u%val) )

end subroutine braid_SpatialNorm_F90


subroutine braid_Access_F90(app, u, astatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8) :: astatus
   type(my_vector)  :: u
   type(my_app)     :: app
   
   ! Other declarations
   integer          :: iter, level, done, ierr, step, numprocs, rank, out_unit
   double precision :: t
   character(len=25):: fname = "ex-01b-f.out"
   character(len=12):: fname_short = "ex-01b-f.out"
   character(len=4) :: step_string
   character(len=3) :: rank_string
   character(len=1) :: dot = "."
   character(len=21) :: val_string
   
   call braid_access_status_get_tild_f90(astatus, t, iter, level, done)
   call braid_access_status_get_tindex_f90(astatus, step)

   ! Print my rank 
   call mpi_comm_size(app%comm, numprocs, ierr)
   call mpi_comm_rank(app%comm, rank, ierr)

   ! Print info to screan
   !print *, "   access print vector", u%val, '  t = ',  t
   !print *, "u->val          = ", u%val
   
   ! Write the scalar solution at this time to file
   write(step_string, "(I4)")  step
   call replace(step_string, ' ', '0', 4)
   write(rank_string, "(I3)")  rank
   call replace(rank_string, ' ', '0', 3)
   write(val_string, "(E21.15)")  u%val
   fname = fname_short // dot // step_string // dot // rank_string

   out_unit = 11
   open (unit=out_unit, file=fname, action="write", status="replace")
   write (out_unit,*) val_string
   close (out_unit)

end subroutine braid_Access_F90


subroutine braid_Step_F90(app, ustop, fstop, fnotzero, u, pstatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8) :: pstatus
   type(my_vector)  :: ustop
   type(my_vector)  :: fstop
   integer          :: fnotzero
   type(my_vector)  :: u
   type(my_app)     :: app

   ! Other declarations
   double precision tstart, tstop, dt

   ! query the status structure for tstart and tstop 
   call braid_step_status_get_tstart_tstop_f90(pstatus, tstart, tstop)
   dt = tstop - tstart

   ! On the finest grid, each value is half the previous value
   u%val = (0.5**dt)*(u%val)

   if (fnotzero .eq. 1) then
      ! Nonzero rhs
      u%val = u%val + fstop%val
   end if

   ! no refinement
   call braid_step_status_set_rfactor_f90(pstatus, 1)

end subroutine braid_Step_F90


subroutine braid_Residual_F90(app, ustop, r, pstatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8) :: pstatus
   type(my_vector)  :: ustop
   type(my_vector)  :: r
   type(my_app)     :: app

   ! Other declarations
   double precision tstart, tstop, dt

   ! query the status structure for tstart and tstop 
   call braid_step_status_get_tstart_tstop_f90(pstatus, tstart, tstop)
   dt = tstop - tstart

   ! On the finest grid, each value is half the previous value
   r%val = (ustop%val) - (0.5**dt)*(r%val)

   ! no refinement
   call braid_step_status_set_rfactor_f90(pstatus, 1)

end subroutine braid_Residual_F90


subroutine braid_BufSize_F90(app, size_ptr, bstatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8)         :: bstatus
   type(my_app)             :: app
   
   ! Other declarations
   integer size_ptr
   size_ptr = (app%bufsize)*sizeof_double

end subroutine braid_BufSize_F90


subroutine braid_BufPack_F90(app, u, buffer, bstatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8)         :: bstatus
   type(my_vector)          :: u
   type(my_app)             :: app
   
   ! Other declarations
   double precision, dimension(app%bufsize) :: buffer
   
   ! Pack buffer
   buffer(1) = u%val
   call braid_buffer_status_set_size_f90(bstatus, app%bufsize*sizeof_double)

end subroutine braid_BufPack_F90


subroutine braid_BufUnPack_F90(app, buffer, u_ptr, bstatus)
   
   ! Braid types
   use braid_types
   implicit none
   integer (kind=8)         :: bstatus
   type(my_vector), pointer :: u_ptr
   type(my_app)             :: app
   
   ! Other declarations
   double precision, dimension(app%bufsize) :: buffer
   
   ! UnPack buffer
   allocate(u_ptr)
   u_ptr%val = buffer(1)

end subroutine braid_BufUnPack_F90


program ex01_f90
   
   ! Import braid vector and app module
   use braid_types
   
   !MPI module for newer f90 compilers. For f77 code, comment in the "include 'mpif.h'" option after "implicit none"  
   !use of "include 'mpif.h'" with f90 compilers may result in "invalid communicator"-type errors
   !use mpi 

   implicit none
   type(my_app), pointer :: app

   ! MPI module for older (f77) compilers; for f90 code "use mpi" at beginning of "program ex01_f90"
   include 'mpif.h'
   
   ! Declare variables
   double precision t, tol
   integer ierr, rc, max_levels, nrelax, nrelax0, cfactor, cfactor0
   integer max_iter, fmg, wrapper_tests, print_help, i, numarg
   integer min_coarse, print_level, access_level, nfmg_Vcyc, res
   character (len = 255) arg
   
   ! Seed random number generator
   call random_seed()

   ! Initialize mpi
   call mpi_init(ierr)
   if (ierr .ne. mpi_success) then
      print *,'Error starting mpi program. Terminating.'
      call mpi_abort(mpi_comm_world, rc, ierr)
   end if

   ! Set default values
   allocate(app)
   app%tstart    = 0.0
   app%ntime     = 10
   app%tstop     = app%tstart + app%ntime
   app%mydt      = 0
   app%bufsize   = 1
   app%comm      = mpi_comm_world
   max_levels    = 2
   nrelax        = 1
   nrelax0       = -1
   tol           = 1.0e-06
   cfactor       = 2
   cfactor0      = -1
   max_iter      = 100
   fmg           = 0
   nfmg_Vcyc     = 1
   res           = 0
   print_help    = 0
   wrapper_tests = 0
   min_coarse    = 3
   print_level   = 1
   access_level  = 1
   
   call mpi_comm_size(app%comm, app%pt,  ierr)
   call mpi_comm_rank(app%comm, app%rank, ierr)

   ! Parse command line
   ! GNU, iFort and PGI support the iargc( ) and getarg( ) functions
   numarg = iargc( )
   i = 0
   do while (i <= numarg)
      call getarg ( i, arg)
      if (arg == '-nt') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) app%ntime
      else if (arg == '-ml') then
         i = i+1; call getarg ( i, arg); i = i+1 
         read(arg,*) max_levels
      else if (arg == '-mc') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) min_coarse
      else if (arg == '-nu') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) nrelax
      else if (arg == '-nu0') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) nrelax0
      else if (arg == '-tol') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) tol
      else if (arg == '-cf') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) cfactor
      else if (arg == '-cf0') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) cfactor0
      else if (arg == '-mi') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) max_iter
      else if (arg == '-fmg') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) nfmg_Vcyc
         fmg = 1
      else if (arg == '-res') then
         i = i+1; call getarg ( i, arg); i = i+1
         res = 1
      else if (arg == '-print_level') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) print_level
      else if (arg == '-access_level') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) access_level
      else if (arg == '-tg') then
         i = i+1; call getarg ( i, arg); i = i+1
         read(arg,*) app%mydt
      else if (arg == '-wrapper_tests') then
         i = i+1
         wrapper_tests = 1
      else if (arg == '-help') then
         i = i+1
         print *, " " 
         print *, "  -nt  <n>            : number of time steps (default: 32)"
         print *, "  -ml  <max_levels>   : set max levels"
         print *, "  -mc  <min_coarse>   : set min possible coarse level size (default: 3)"
         print *, "  -nu  <nrelax>       : set num F-C relaxations"
         print *, "  -nu0 <nrelax>       : set num F-C relaxations on level 0"
         print *, "  -tol <tol>          : set stopping tolerance"
         print *, "  -cf  <cfactor>      : set coarsening factor"
         print *, "  -cf0  <cfactor>     : set coarsening factor for level 0 "
         print *, "  -mi  <max_iter>     : set max iterations"
         print *, "  -fmg <nfmg_Vcyc>    : use FMG cycling, nfmg_Vcyc  V-cycles at each fmg level"
         print *, "  -res                : use user-defined residual computation function"
         print *, "                        must compile with correct option in braid.h "
         print *, "  -tg <mydt>          : use user-specified time grid as global fine time grid, options are"
         print *, "                        1 - uniform time grid"
         print *, "                        2 - nonuniform time grid, dt*0.5 for n = 1, ..., nt/2; dt*1.5 for n = nt/2+1, ..., nt"
         print *, "  -print_level <l>    : sets the print_level (default: 1) "
         print *, "                        0 - no output to standard out "
         print *, "                        1 - Basic convergence information and hierarchy statistics"
         print *, "                        2 - Debug level output"
         print *, "  -access_level <l>   : sets the access_level (default: 1)"
         print *, "                        0 - never call access"
         print *, "                        1 - call access only after completion"
         print *, "                        2 - call access every iteration and level"
         print *, "  -wrapper_tests      : Only run the XBraid wrapper tests"
         print *, " "
         print_help = 1
         exit
      else
         i = i + 1
      endif
   enddo
   
   if (print_help == 0) then
      if (wrapper_tests == 1) then
         ! Call Braid Test Routines
         t = 0.1
         call braid_test_all_f90(app, mpi_comm_world, t, t, t)
      else

         ! get rank for local communicators
         call mpi_comm_rank(app%comm, app%rank, ierr)

         ! Initialize braid
         call braid_init_f90(app%comm, app%comm, app%tstart, app%tstop, &
              app%ntime, app, braid_core)

         ! Set the various command line options for braid
         call braid_set_max_levels_f90(braid_core, max_levels)
         call braid_set_nrelax_f90(braid_core, -1, nrelax)
         if (nrelax0 > -1) then
            call braid_set_nrelax_f90(braid_core,  0, nrelax0)
         endif
         call braid_set_abs_tol_f90(braid_core, tol)
         call braid_set_cfactor_f90(braid_core, -1, cfactor)
         if (cfactor0 > -1) then
            call braid_set_cfactor_f90(braid_core,  0, cfactor0)
         endif
         call braid_set_max_iter_f90(braid_core, max_iter)
         if (fmg == 1) then
            call braid_set_fmg_f90(braid_core)
            call braid_set_nfmg_vcyc_f90(braid_core, nfmg_Vcyc)
         endif
         call braid_set_min_coarse_f90(braid_core, min_coarse)
         call braid_set_print_level_f90( braid_core, print_level)
         call braid_set_access_level_f90( braid_core, access_level)

         ! Define the initial fine time grid using a user supplied function
         if (app%mydt > 0) then
            call init_timesteps(app)
            call braid_set_timegrid_f90(braid_core)
         endif
         
         ! Run braid and clean up
         call braid_drive_f90(braid_core)
         call braid_destroy_f90(braid_core)
         
      endif
   endif

   ! Clean up
   if (allocated(app%dt)) deallocate(app%dt)
   deallocate(app)
   call mpi_finalize(ierr)

end program ex01_f90

