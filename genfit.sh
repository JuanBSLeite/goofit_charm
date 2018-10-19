#!/bin/sh

for i in $(seq 0 100); do
	./D2PPP -i $i gen
done

for i in $(seq 0 100); do
	./D2PPP -i $i fit	
done

./D2PPP gfplot -n 100