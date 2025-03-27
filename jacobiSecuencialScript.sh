#!/bin/bash

output_file="resultadoJacobiSecuencial_O3.csv"

# Comprueba si el archivo ya existe y, si no, añade los encabezados
if [ ! -f "$output_file" ]; then
  echo "Longitud,Tiempo" > "$output_file"
fi

for i in 50 100 1000 3000 5000; do
  echo "Ejecutando script para arreglo de longitud: $i"

  for j in {1..10}; do
    ./jacobiSecuencialExe_O3 $i $i >> "$output_file"
  done
    
  echo "" >> "$output_file"
done

echo "Pruebas finalizadas."
