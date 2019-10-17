#!/bin/bash

echo "Script to run simulation as a function of T"
echo "Removing old files..."

rm ene.dat
rm mag.dat
rm heat.dat
rm chi.dat

TEMPS=($(seq 0.5 0.05 2.05))

for temp in "${TEMPS[@]}"
do
    echo "Simulating T=$temp K ..."
    gsed -i "1s/.*/$temp/g" input.dat  #gsed for macos, on gnu use 'sed'
    ./clean.sh
    ./Monte_Carlo_ISING_1D.exe > /dev/null
    tail -1 output.ene.0 | awk -v var="$temp" '{print var, $3, $4}' >> ene.dat
    tail -1 output.mag.0 | awk -v var="$temp" '{print var, $3, $4}' >> mag.dat
    tail -1 output.heat.0  | awk -v var="$temp" '{print var, $3, $4}' >> heat.dat
    tail -1 output.chi.0  | awk -v var="$temp" '{print var, $3, $4}' >> chi.dat
done

echo "All done!"
