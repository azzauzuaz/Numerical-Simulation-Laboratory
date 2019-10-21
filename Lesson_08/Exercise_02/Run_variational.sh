#!/bin/bash

MU_MIN=0.7
MU_MAX=1.0
MU_STEP=0.01

SIGMA_MIN=0.5
SIGMA_MAX=0.8
SIGMA_STEP=0.01

echo "Script to run Variational Monte Carlo"

echo "Generating parameters..."

for mu in `seq $MU_MIN $MU_STEP $MU_MAX`
do
  for sigma in `seq $SIGMA_MIN $SIGMA_STEP $SIGMA_MAX`
  do
    echo "$mu $sigma" >> parameters
  done
done

if [ -f variational.dat ]; then
    rm variational.dat
fi

echo "Running simulation..."

while read mu sigma
do
    echo "Simulating mu=$mu, sigma=$sigma..."
    gsed -i "5s/.*/$mu/g" input.dat  #gsed for macos, on gnu use 'sed'
    gsed -i "6s/.*/$sigma/g" input.dat
    rm hamiltonian.dat
    ./metropolis.exe > /dev/null
    tail -1 hamiltonian.dat | awk -v mu="$mu" -v sigma="$sigma" '{print mu, sigma, $2, $3}' >> variational.dat

done < parameters

rm parameters

export LC_ALL=C #mandatory on macos to make sure awk and sort don't play strange jokes (shoud be harmless on gnu)

RESULTS=($(awk '{print $3, $1, $2}' variational.dat | sort -n | head -1))

echo "Optimal parameters are: mu=${RESULTS[1]}, sigma=${RESULTS[2]}"

echo "Results saved in 'variational.dat'..."

echo "All done!"
