#!/bin/sh
START_TIME=$SECONDS
#for i in $(seq 0 99); do
#	./D2PPP -i $i gen -e 100000
#done

for i in $(seq 0 99); do
	./D2PPP -i $i fit	
done

./D2PPP gfplot -n 100 -v 10

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   


 
