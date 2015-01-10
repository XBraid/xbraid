srun -N 2 -n 24 ./drive-05 -pgrid 1 1 24 -nt 2048 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 4 -scoarsen 2 -scoarsenCFL 2 10.0 5.0

srun -N 2 -n 16  ./drive-05 -pgrid 1 1 16 -nt 2048 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 4 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3

# Re-runs problem from figure 4.4 in paper, where cross over was 16 processors in time
# --> Is residual tolerance too tight?
# --> Is 257x257 on  each proc just paging stuff out of mem?

# SHOULD THE FIRST CFL NUMBER BE LIKE 100.0 ??  REMEMBER, NO SPATIAL COARSENING

srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 4 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3
srun -N 1 -n 4 ./drive-05 -pgrid 2 2 1 -nt 16385 -ml 1 -nx 129 129 
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 8.226448e+01,  wall time = 4.41e+00
  Braid:  || r_1 || = 8.672833e+00,  wall time = 8.77e+00,  conv. factor = 1.05e-01
  Braid:  || r_2 || = 8.315300e-01,  wall time = 1.31e+01,  conv. factor = 9.59e-02
  Braid:  || r_3 || = 8.229161e-02,  wall time = 1.75e+01,  conv. factor = 9.90e-02
  Braid:  || r_4 || = 8.224131e-03,  wall time = 2.18e+01,  conv. factor = 9.99e-02
  Braid:  || r_5 || = 8.289613e-04,  wall time = 2.62e+01,  conv. factor = 1.01e-01
  Braid:  || r_6 || = 8.427038e-05,  wall time = 3.05e+01,  conv. factor = 1.02e-01
  Braid:  || r_7 || = 8.613865e-06,  wall time = 3.49e+01,  conv. factor = 1.02e-01
  Braid:  || r_8 || = 8.826145e-07,  wall time = 3.93e+01,  conv. factor = 1.02e-01
  my_Access() called, iter= 9, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 5
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 9
  residual norm        = 8.826145e-07
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256        4        1
      2        64        4        1
      3        16        4        1
      4         4  

  wall time = 41.460747

  runtime: 41.46149s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   4.91e-02    4.91e-02    5.78e-03    4.80e+00
  2   |  impl,  3   9.82e-02    9.82e-02    2.31e-02    4.80e+00
  3   |  impl,  3   1.96e-01    1.96e-01    9.25e-02    4.80e+00
  4   |  impl,  3   3.93e-01    3.93e-01    3.70e-01    4.80e+00


srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 16 -cf 4 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 1.268814e+02,  wall time = 5.03e+00
  Braid:  || r_1 || = 1.480623e+01,  wall time = 9.93e+00,  conv. factor = 1.17e-01
  Braid:  || r_2 || = 1.508270e+00,  wall time = 1.48e+01,  conv. factor = 1.02e-01
  Braid:  || r_3 || = 1.587399e-01,  wall time = 1.97e+01,  conv. factor = 1.05e-01
  Braid:  || r_4 || = 1.687665e-02,  wall time = 2.46e+01,  conv. factor = 1.06e-01
  Braid:  || r_5 || = 1.810190e-03,  wall time = 2.96e+01,  conv. factor = 1.07e-01
  Braid:  || r_6 || = 1.954473e-04,  wall time = 3.45e+01,  conv. factor = 1.08e-01
  Braid:  || r_7 || = 2.120376e-05,  wall time = 3.94e+01,  conv. factor = 1.08e-01
  Braid:  || r_8 || = 2.308362e-06,  wall time = 4.43e+01,  conv. factor = 1.09e-01
  my_Access() called, iter= 9, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 6
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 9
  residual norm        = 2.308362e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       16        1
      1      1024        4        1
      2       256        4        1
      3        64        4        1
      4        16        4        1
      5         4  

  wall time = 46.270349

  runtime: 46.27041s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   2.45e-02    2.45e-02    1.45e-03    4.80e+00
  2   |  impl,  3   4.91e-02    4.91e-02    5.78e-03    4.80e+00
  3   |  impl,  3   9.82e-02    9.82e-02    2.31e-02    4.80e+00
  4   |  impl,  3   1.96e-01    1.96e-01    9.25e-02    4.80e+00
  5   |  impl,  3   3.93e-01    3.93e-01    3.70e-01    4.80e+00



srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 8.226448e+01,  wall time = 4.41e+00
  Braid:  || r_1 || = 7.008633e+00,  wall time = 8.75e+00,  conv. factor = 8.52e-02
  Braid:  || r_2 || = 5.303350e-01,  wall time = 1.31e+01,  conv. factor = 7.57e-02
  Braid:  || r_3 || = 4.496194e-02,  wall time = 1.74e+01,  conv. factor = 8.48e-02
  Braid:  || r_4 || = 3.975407e-03,  wall time = 2.18e+01,  conv. factor = 8.84e-02
  Braid:  || r_5 || = 3.592740e-04,  wall time = 2.61e+01,  conv. factor = 9.04e-02
  Braid:  || r_6 || = 3.282043e-05,  wall time = 3.05e+01,  conv. factor = 9.14e-02
  Braid:  || r_7 || = 3.017700e-06,  wall time = 3.48e+01,  conv. factor = 9.19e-02
  my_Access() called, iter= 8, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 8
  residual norm        = 3.017700e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 36.917448

  runtime: 36.91751s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   4.91e-02    4.91e-02    5.78e-03    4.80e+00
  2   |  impl,  3   1.96e-01    1.96e-01    9.25e-02    4.80e+00



srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 1 0.4 -iter 50 3
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 8.226448e+01,  wall time = 4.35e+00
  Braid:  || r_1 || = 1.282950e+01,  wall time = 8.67e+00,  conv. factor = 1.56e-01
  Braid:  || r_2 || = 1.719563e+00,  wall time = 1.30e+01,  conv. factor = 1.34e-01
  Braid:  || r_3 || = 3.008760e-01,  wall time = 1.73e+01,  conv. factor = 1.75e-01
  Braid:  || r_4 || = 6.360293e-02,  wall time = 2.17e+01,  conv. factor = 2.11e-01
  Braid:  || r_5 || = 1.475864e-02,  wall time = 2.60e+01,  conv. factor = 2.32e-01
  Braid:  || r_6 || = 3.624865e-03,  wall time = 3.03e+01,  conv. factor = 2.46e-01
  Braid:  || r_7 || = 9.250882e-04,  wall time = 3.46e+01,  conv. factor = 2.55e-01
  Braid:  || r_8 || = 2.424032e-04,  wall time = 3.89e+01,  conv. factor = 2.62e-01
  Braid:  || r_9 || = 6.471243e-05,  wall time = 4.33e+01,  conv. factor = 2.67e-01
  Braid:  || r_10 || = 1.751453e-05,  wall time = 4.76e+01,  conv. factor = 2.71e-01
  Braid:  || r_11 || = 4.790696e-06,  wall time = 5.19e+01,  conv. factor = 2.74e-01
  Braid:  || r_12 || = 1.321414e-06,  wall time = 5.63e+01,  conv. factor = 2.76e-01
  my_Access() called, iter= 13, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 13
  residual norm        = 1.321414e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 58.377069

  runtime: 58.37712s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   1.96e-01    1.96e-01    5.78e-03    3.00e-01
  2   |  impl,  3   7.85e-01    7.85e-01    9.25e-02    3.00e-01



srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 16 -iter 50 3

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 8.226448e+01,  wall time = 4.46e+00
  Braid:  || r_1 || = 8.341324e+00,  wall time = 8.91e+00,  conv. factor = 1.01e-01
  Braid:  || r_2 || = 7.123532e-01,  wall time = 1.34e+01,  conv. factor = 8.54e-02
  Braid:  || r_3 || = 6.834217e-02,  wall time = 1.78e+01,  conv. factor = 9.59e-02
  Braid:  || r_4 || = 6.883841e-03,  wall time = 2.23e+01,  conv. factor = 1.01e-01
  Braid:  || r_5 || = 7.138010e-04,  wall time = 2.68e+01,  conv. factor = 1.04e-01
  Braid:  || r_6 || = 7.514857e-05,  wall time = 3.12e+01,  conv. factor = 1.05e-01
  Braid:  || r_7 || = 7.975215e-06,  wall time = 3.57e+01,  conv. factor = 1.06e-01
  Braid:  || r_8 || = 8.505708e-07,  wall time = 4.01e+01,  conv. factor = 1.07e-01
  my_Access() called, iter= 9, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 9
  residual norm        = 8.505708e-07
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 42.270841

  runtime: 42.27090s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   2.45e-02    2.45e-02    5.78e-03    1.92e+01
  2   |  impl,  3   2.45e-02    2.45e-02    9.25e-02    3.07e+02



srun -N 1 -n 1 ./drive-05 -pgrid 1 1 1 -nt 16385 -ml 1 -nx 129 129 
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  my_Access() called, iter= 0, level= 0



  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 1
  min coarse           = 3
  number of levels     = 1
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 100
  iterations           = 0
  residual norm        = -1.000000e+00
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385  

  wall time = 32.979732

  runtime: 32.97978s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  50

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  3   2.45e-02    2.45e-02    9.04e-05    3.00e-01








srun -N 1 -n 1 ./drive-05 -pgrid 1 1 1 -nt 16385 -ml 1 -nx 257 257 

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  my_Access() called, iter= 0, level= 0

  start time = 0.000000e+00
  stop time  = 3.701328e-01
  time steps = 16385

  max number of levels = 1
  min coarse           = 3
  number of levels     = 1
  stopping tolerance   = 1.714488e-05
  relative tolerance?  = 0
  max iterations       = 100
  iterations           = 0
  residual norm        = -1.000000e+00
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385  

  wall time = 88.830170

  runtime: 88.83023s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  50

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 257 x 257

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  2   1.23e-02    1.23e-02    2.26e-05    3.00e-01



srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 257 257 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 


  Braid:  || r_0 || = 1.308890e+02,  wall time = 1.83e+01
  Braid:  || r_1 || = 1.033353e+01,  wall time = 3.64e+01,  conv. factor = 7.89e-02
  Braid:  || r_2 || = 7.755199e-01,  wall time = 4.96e+01,  conv. factor = 7.50e-02
  Braid:  || r_3 || = 6.543707e-02,  wall time = 6.27e+01,  conv. factor = 8.44e-02
  Braid:  || r_4 || = 5.804191e-03,  wall time = 7.58e+01,  conv. factor = 8.87e-02
  Braid:  || r_5 || = 5.264194e-04,  wall time = 8.89e+01,  conv. factor = 9.07e-02
  Braid:  || r_6 || = 4.821170e-05,  wall time = 1.02e+02,  conv. factor = 9.16e-02
  Braid:  || r_7 || = 4.442205e-06,  wall time = 1.15e+02,  conv. factor = 9.21e-02
  my_Access() called, iter= 8, level= 0

  start time = 0.000000e+00
  stop time  = 3.701328e-01
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 1.714488e-05
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 8
  residual norm        = 4.442205e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 121.631692

  runtime: 121.63175s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 257 x 257

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   1.23e-02    1.23e-02    2.26e-05    3.00e-01
  1   |  impl,  3   2.45e-02    2.45e-02    1.45e-03    4.80e+00
  2   |  impl,  3   9.82e-02    9.82e-02    2.31e-02    4.80e+00



srun -N 2 -n 16 ./drive-05 -pgrid 1 1 16 -nt 16385 -mi 75 -ml 15 -nx 257 257 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 100.0 5.0 -iter 50 3
  Braid:  || r_0 || = 1.308890e+02,  wall time = 1.88e+01
  Braid:  || r_1 || = 1.136976e+01,  wall time = 3.74e+01,  conv. factor = 8.69e-02
  Braid:  || r_2 || = 9.994712e-01,  wall time = 5.09e+01,  conv. factor = 8.79e-02
  Braid:  || r_3 || = 9.864059e-02,  wall time = 6.45e+01,  conv. factor = 9.87e-02
  Braid:  || r_4 || = 1.015311e-02,  wall time = 7.80e+01,  conv. factor = 1.03e-01
  Braid:  || r_5 || = 1.063635e-03,  wall time = 9.15e+01,  conv. factor = 1.05e-01
  Braid:  || r_6 || = 1.124101e-04,  wall time = 1.05e+02,  conv. factor = 1.06e-01
  Braid:  || r_7 || = 1.194927e-05,  wall time = 1.19e+02,  conv. factor = 1.06e-01
  my_Access() called, iter= 8, level= 0

  start time = 0.000000e+00
  stop time  = 3.701328e-01
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 1.714488e-05
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 8
  residual norm        = 1.194927e-05
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 124.996758

  runtime: 124.99737s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 257 x 257

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   1.23e-02    1.23e-02    2.26e-05    3.00e-01
  1   |  impl,  3   1.23e-02    1.23e-02    1.45e-03    1.92e+01
  2   |  impl,  3   9.82e-02    9.82e-02    2.31e-02    4.80e+00











srun -N 1 -n 4 ./drive-05 -pgrid 2 2 1 -nt 16385 -ml 1 -nx 129 129 

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  my_Access() called, iter= 0, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 1
  min coarse           = 3
  number of levels     = 1
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 100
  iterations           = 0
  residual norm        = -1.000000e+00
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385  

  wall time = 15.758898

  runtime: 15.75894s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  50

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  3   2.45e-02    2.45e-02    9.04e-05    3.00e-01




srun -N 2 -n 32 ./drive-05 -pgrid 2 2 8 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3
  Begin simulation 

  Braid:  || r_0 || = 8.314331e+01,  wall time = 4.16e+00
  Braid:  || r_1 || = 7.080821e+00,  wall time = 8.39e+00,  conv. factor = 8.52e-02

  Braid:  || r_2 || = 5.398875e-01,  wall time = 1.25e+01,  conv. factor = 7.62e-02
  Braid:  || r_3 || = 4.601864e-02,  wall time = 1.67e+01,  conv. factor = 8.52e-02
  Braid:  || r_4 || = 4.084594e-03,  wall time = 2.09e+01,  conv. factor = 8.88e-02
  Braid:  || r_5 || = 3.701409e-04,  wall time = 2.50e+01,  conv. factor = 9.06e-02
  Braid:  || r_6 || = 3.388132e-05,  wall time = 2.92e+01,  conv. factor = 9.15e-02
  Braid:  || r_7 || = 3.120094e-06,  wall time = 3.33e+01,  conv. factor = 9.21e-02
  my_Access() called, iter= 8, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 8
  residual norm        = 3.120094e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 35.292194

  runtime: 35.29292s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   4.91e-02    4.91e-02    5.78e-03    4.80e+00
  2   |  impl,  3   1.96e-01    1.96e-01    9.25e-02    4.80e+00






srun -N 2 -n 32 ./drive-05 -pgrid 2 1 16 -nt 16385 -mi 75 -ml 15 -nx 129 129 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3
-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 8.278441e+01,  wall time = 3.08e+00
  Braid:  || r_1 || = 7.033749e+00,  wall time = 6.14e+00,  conv. factor = 8.50e-02
  Braid:  || r_2 || = 5.324835e-01,  wall time = 9.20e+00,  conv. factor = 7.57e-02
  Braid:  || r_3 || = 4.513083e-02,  wall time = 1.22e+01,  conv. factor = 8.48e-02
  Braid:  || r_4 || = 3.987714e-03,  wall time = 1.53e+01,  conv. factor = 8.84e-02


  Braid:  || r_5 || = 3.600969e-04,  wall time = 1.84e+01,  conv. factor = 9.03e-02
  Braid:  || r_6 || = 3.286942e-05,  wall time = 2.14e+01,  conv. factor = 9.13e-02
  Braid:  || r_7 || = 3.019984e-06,  wall time = 2.45e+01,  conv. factor = 9.19e-02
  my_Access() called, iter= 8, level= 0

  start time = 0.000000e+00
  stop time  = 1.480531e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 4.286221e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 8
  residual norm        = 3.019984e-06
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 26.006126

  runtime: 26.00694s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 129 x 129

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   2.45e-02    2.45e-02    9.04e-05    3.00e-01
  1   |  impl,  3   4.91e-02    4.91e-02    5.78e-03    4.80e+00
  2   |  impl,  3   1.96e-01    1.96e-01    9.25e-02    4.80e+00









cab292@schroder:srun -N 1 -n 4 ./drive-05 -pgrid 2 2 1 -nt 16385 -ml 1 -nx 65 65 

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  my_Access() called, iter= 0, level= 0

  start time = 0.000000e+00
  stop time  = 5.922124e+00
  time steps = 16385

  max number of levels = 1
  min coarse           = 3
  number of levels     = 1
  stopping tolerance   = 1.071555e-06
  relative tolerance?  = 0
  max iterations       = 100
  iterations           = 0
  residual norm        = -1.000000e+00
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385  

  wall time = 8.344230

  runtime: 8.34427s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  50

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 65 x 65

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  3   4.91e-02    4.91e-02    3.61e-04    3.00e-01




cab292@schroder:srun -N 2 -n 32 ./drive-05 -pgrid 1 1 32 -nt 16385 -mi 75 -ml 15 -nx 65 65 -cf0 64 -cf 16 -scoarsen 2 -scoarsenCFL 2 10.0 5.0 -iter 50 3

-----------------------------------------------------------------
-----------------------------------------------------------------

  Begin simulation 

  Braid:  || r_0 || = 5.518179e+01,  wall time = 6.92e-01
  Braid:  || r_1 || = 5.264023e+00,  wall time = 1.39e+00,  conv. factor = 9.54e-02
  Braid:  || r_2 || = 3.837009e-01,  wall time = 2.08e+00,  conv. factor = 7.29e-02
  Braid:  || r_3 || = 3.179849e-02,  wall time = 2.77e+00,  conv. factor = 8.29e-02
  Braid:  || r_4 || = 2.821196e-03,  wall time = 3.47e+00,  conv. factor = 8.87e-02
  Braid:  || r_5 || = 2.585139e-04,  wall time = 4.16e+00,  conv. factor = 9.16e-02
  Braid:  || r_6 || = 2.400443e-05,  wall time = 4.86e+00,  conv. factor = 9.29e-02
  Braid:  || r_7 || = 2.245451e-06,  wall time = 5.55e+00,  conv. factor = 9.35e-02
  Braid:  || r_8 || = 2.109826e-07,  wall time = 6.24e+00,  conv. factor = 9.40e-02
  my_Access() called, iter= 9, level= 0

  start time = 0.000000e+00
  stop time  = 5.922124e+00
  time steps = 16385

  max number of levels = 15
  min coarse           = 3
  number of levels     = 3
  stopping tolerance   = 1.071555e-06
  relative tolerance?  = 0
  max iterations       = 75
  iterations           = 9
  residual norm        = 2.109826e-07
                        --> 2-norm TemporalNorm 

  level   time-pts   cfactor   nrelax
      0     16385       64        1
      1       256       16        1
      2        16  

  wall time = 6.584298

  runtime: 6.58507s


-----------------------------------------------------------------
-----------------------------------------------------------------

 Implicit time stepping solve parameters

   Fine-level loose stopping tol  :  1.00e-09    (while ||r|| is large)
   Fine-level tight stopping tol  :  1.00e-09    (while ||r|| is small)
   Coarse-level stopping tol      :  1.00e-09    (for all ||r||) 
 
   Fine-level max iter            :  50
   Coarse-level max iter          :  3

-----------------------------------------------------------------
-----------------------------------------------------------------

 Per level diagnostic information 

 Scheme is either explicit or implicit.  If implicit then the table
 entry is 'impl, %d' where the integer represents the maximum number
 of AMG iterations used by the implicit solver on that level.  'expl'
 stands for explicit. 

 Fine level spatial problem size  : 65 x 65

level   scheme         dx          dy          dt          cfl
-----------------------------------------------------------------
  0   |  impl,  4   4.91e-02    4.91e-02    3.61e-04    3.00e-01
  1   |  impl,  3   9.82e-02    9.82e-02    2.31e-02    4.80e+00
  2   |  impl,  3   3.93e-01    3.93e-01    3.70e-01    4.80e+00

 
