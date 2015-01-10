# Coarsen more aggressively in space (CFL decreases by a factor of 2)
mpirun -np 4 drive-05 -pgrid 1 1 4 -nt 4096 -nx 65 65 -ml 2 -scoarsen 1 -cf0 2 -nu0 1 -cf 2 -cfl 0.48 -scoarsenCFL 0.25 -expl -mi 30


# Coarsen more aggressively in space (CFL decreases by a factor of 8)
mpirun -np 4 drive-05 -pgrid 1 1 4 -nt 4096 -nx 65 65 -ml 2 -scoarsen 1 -cf0 2 -nu0 1 -cf 2 -cfl 0.48 -scoarsenCFL 0.1 -expl -mi 30


# rm drive-05.out.iter* drive-05_mesh.* drive-05_sol* drive-05_err*

