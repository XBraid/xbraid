fargs="-tf 4 -nx 128 -nt 2048 -theta -out -Delta -rank 18 -Deltalvl 0 -tol 1e-4"

echo "mpirun drive-ks ${fargs} -ml 1"
mpirun drive-ks ${fargs} -ml 1 -cglv

echo "mv drive-ks.out drive-ks-seq.out"
echo "mv drive-ks-lv.out drive-ks-lv-seq.out"
mv drive-ks.out drive-ks-seq.out
mv drive-ks-lv.out drive-ks-lv-seq.out

echo "mpirun drive-ks ${fargs}"
mpirun drive-ks ${fargs} -ml 8
