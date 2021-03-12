#!/bin/bash
for i in $(seq 2 1000); do
    ./D2PPP_rhoStudies --gpu-dev=1 -q fit -n rhoStudies/Fit_$i -f ../../../dados/DsPPP_Reduced_93.root
    echo "Fit: $i"
done
