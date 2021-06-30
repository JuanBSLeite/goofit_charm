#!/bin/bash
M=0
for i in $(seq 1.416 0.025 1.516); do
	for j in $(seq 0.28 0.06 0.52); do
		for k in $(seq 1.68 0.02 1.76); do
			for l in $(seq 0.2 0.1 0.4); do
				./D2PPP_model1_scan --gpu-dev=1 -q fit -n model1_resonance_scan/Fit_$M --r1-Mass $i --r1-Width $j --r2-Mass $k --r2-Width $l -f ../../../8-FitSamples/DsPPP_SignalRegion_AllSelectionCuts_Splot_FReweighted_MeanMVA_95.root -s 1
				let M+=1
				echo "M = $M"
			done
		done
	done
done
