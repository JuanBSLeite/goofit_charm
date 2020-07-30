#!/bin/sh
START_TIME=$SECONDS

for i in $(seq 0 9); do
	./D2PPP -q genfit -n $i 
done

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   


 
