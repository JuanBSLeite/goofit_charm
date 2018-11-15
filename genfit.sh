#!/bin/sh

for i in $(seq 0 2); do
	./D2PPP -i $i gen -e 250000
done

for i in $(seq 0 2); do
	./D2PPP -i $i fit	
done

./D2PPP gfplot -n 3