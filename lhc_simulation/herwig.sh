#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DEBUG=0

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))

if [ $JOB_TYPE -eq 0 ]; then
    HFILE=LHC-dijet
elif [ $JOB_TYPE -eq 1 ]; then
    HFILE=LHC-dijet
elif [ $JOB_TYPE -eq 2 ]; then
    HFILE=LHC-gammajet
elif [ $JOB_TYPE -eq 3 ]; then
    HFILE=LHC-Zjet
else
    HFILE=LHC-dijet
fi

cd herwig
Herwig++ read ${HFILE}.in
wait $!

pidArr=()
for (( i=1; i<=$NUM_PROC; i++ ))
do
    Herwig++ run -N$EVT_PER_RUN -s$((100*$i)) -t++${NUM_PROC}_${i}_${JOB_TYPE} ${HFILE}.run &
    pidArr+=($!)
    pidArr+=" "
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done

MERGE=$(ls -rt | grep root | tail -n $(($NUM_PROC+1)) | head -n 1)
TEMPORARY=$(ls -rt | grep root | tail -n $NUM_PROC)

if [ $NUM_PROC -gt 1 ]; then
    hadd -f $MERGE $TEMPORARY

    for tmp in $TEMPORARY
    do
        rm $tmp
    done
else
    MERGE=$TEMPORARY
fi

if [ $DEBUG -eq 0 ]; then
    rm *.out
    rm *.log
    rm *.tex
fi

mv $MERGE ..

exit
