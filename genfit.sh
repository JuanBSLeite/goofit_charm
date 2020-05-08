#!/bin/sh
START_TIME=$SECONDS


#for i in $(seq 0 99); do
#	./D2PPP -q -f "MC_$i" --make-toy=true --save-toy=true 
#done

for i in $(seq 0 99); do
	./D2PPP -q -f "MC_$i" --genfit=true	
done

#./D2PPP gfplot -n 500 -v 4

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   


 
