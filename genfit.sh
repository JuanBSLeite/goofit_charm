#!/bin/sh
START_TIME=$SECONDS
#for i in $(seq 0 300); do
#	./D2PPP -i $i gen -e 250000
#done

for i in $(seq 40 300); do
	./D2PPP -i $i fit	
done

./D2PPP gfplot -n 301 -v 113

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   


 
