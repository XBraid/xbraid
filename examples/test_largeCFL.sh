CFLS=( "0.45" "0.48" "0.49" "0.499" )
GRIDS=( "-nt 256 -nx 17 17" "-nt 1024 -nx 33 33" "-nt 4096 -nx 65 65" "-nt 16384 -nx 129 129" "-nt 65536 -nx 257 257")
CFS=("2" "4" "8" "16" "32" "64" "128" "256" "512" "1024")

for C in "${CFLS[@]}"; do
   for G in "${GRIDS[@]}"; do
      for CF in "${CFS[@]}"; do
         echo ""
         echo ""
         echo "--------------------------------------------------------------------"
         echo "--------------------------------------------------------------------"
         echo ""
         echo ""
         echo "mpirun -np 6 drive-05 -pgrid 1 1 6 $G -ml 2 -scoarsen 1 -cf0 $CF -nu0 1 -cf 2 -cfl $C -scoarsenCFL 0.5 -expl -mi 500 -print_level 1"
         echo ""
         echo ""
         mpirun -np 6 drive-05 -pgrid 1 1 6 $G -ml 2 -scoarsen 1 -cf0 $CF -nu0 1 -cf 2 -cfl $C -scoarsenCFL 0.5 -expl -mi 500 -print_level 1
      done
   done
done
