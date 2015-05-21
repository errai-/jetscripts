#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))

cd pythia6
pidArr=()
for (( i=1; i<=$NUM_PROC; i++ ))
do
    ./pythia6.exe $EVT_PER_RUN $JOB_TYPE $i $NUM_PROC &
    pidArr+=($!)
    pidArr+=" "
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done

MERGE=$(ls -rt | grep root | tail -n $(($NUM_PROC+1)) | head -n 1)
TEMPORARY=$(ls -rt | grep root | tail -n $(($NUM_PROC)))

if [ $NUM_PROC -gt 1 ]; then
    hadd -f $MERGE $TEMPORARY

    for tmp in $TEMPORARY
    do
        rm $tmp
    done
else
    MERGE=$TEMPORARY
fi

mv $MERGE ..

exit
