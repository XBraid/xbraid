This is an example using the Braid-Cython interface defined in braid.pyx (
braid/braid.pyx ). This example uses a higher-level more Python-style syntax
than the other basic Cython example in `examples/ex-01-cython-alt`
 

It solves the same scalar ODE equation 
 
 u' = lambda u, 
 with lambda=-1 and y(0) = 1

as the ex-01 series in the `examples` directory. 

For instructions on running and compiling, see 

   examples/ex-01-cython/ex_01-setup.py

and 

   examples/ex-01-cython/ex_01.pyx

