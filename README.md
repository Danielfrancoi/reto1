# reto1
reto1 HPC jacobi1D

## SECUENCIAL, O1, O2, O3

gcc jacobiSecuencial.c -o jacobiSecuencialExe

gcc -O1 jacobiSecuencial.c -o jacobiSecuencialExe
gcc -O2 jacobiSecuencial.c -o jacobiSecuencialExe
gcc -O3 jacobiSecuencial.c -o jacobiSecuencialExe

ejecutar pruebas

bash ./jacobiSecuencialScript.sh

## HILOS

gcc -o jacobiHilosExe jacobiHilos.c -pthread

ejecutar pruebas

bash ./jacobiHilosScript.sh

## PROCESOS

gcc -g jacobiProcesos.c -o jacobiProcesosExe

ejecutar pruebas

bash ./jacobiProcesosScript.sh
