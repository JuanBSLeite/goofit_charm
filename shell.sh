START_TIME=$SECONDS

PHA_LIST=(0 5 10 15 20) 
AMP_LIST=(1)  

for MAG in ${AMP_LIST[*]}; do
	for PHS in ${PHA_LIST[*]}; do
		./D2PPP -q  --make-toy=true --save-toy=true --phi-1020-mag $MAG --phi-1020-phs $PHS --nEvents=100000000
	done
done

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"   
