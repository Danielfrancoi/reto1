#!/bin/bash

output_file="resultadoJacobiOpenmp.csv"

# Comprueba si el archivo ya existe y, si no, añade los encabezados
if [ ! -f "$output_file" ]; then
    echo "Longitud,Tiempo,Hilos" > "$output_file"
fi

for i in 10 100 1000 10000 100000; do
  echo "Ejecutando script para arreglo de longitud: $i"

  for j in {1..10}; do
    ./jacobiOpenmp $i $i >> "$output_file"
  done

   echo "" >> "$output_file"
done

echo "Pruebas finalizadas."

    