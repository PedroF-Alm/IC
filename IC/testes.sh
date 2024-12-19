#!/bin/bash

# Bateria de testes

rm -r testes
mkdir testes
mkdir testes/tempos
mkdir testes/graficos
mkdir testes/graficos/ap
mkdir testes/graficos/ical

executar () {
    echo $0
    cd $1
    cmake ..
    make
    { time ./exec 4 16x16 > "../../testes/tempos/$2.txt" 2> /dev/zero; } 2>> "../../testes/tempos/$2.txt"
    cd ../exit/4_lcc_16x16/
    python3 plot.py
    mv ap.pdf "../../../testes/graficos/ap/ap_$2.pdf"
    mv ical.pdf "../../../testes/graficos/ical/ical_$2.pdf"
    cd ../../..
}

for i in {1..5}
do
    executar ./doutorado_final_serial/build/ "serial$i"
    executar ./MPI/doutorado_final_mpi/build/ "mpi$i"
    executar ./OPENMP/doutorado_final_openmp/build/ "openmp$i"
    executar ./CUDA/doutorado_final_cuda/build/ "cuda$i"
done