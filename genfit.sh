#!/bin/sh
START_TIME=$SECONDS


#for i in $(seq 0 99); do
#	./D2PPP -q -f "MC_$i" --make-toy=true --save-toy=true 
#done

#for toy genfit
#for i in $(seq 6 99); do
#	./D2PPP_genfit -q -f "MC_$i" --genfit=true	
#done

#multiple soluctions studies

#./D2PPP gfplot -n 500 -v 4
for i in $(seq 0 2); do
	./D2PPP -q -f ../../../dados/DsPPP_92.root --genfit=true	
done

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   


 
