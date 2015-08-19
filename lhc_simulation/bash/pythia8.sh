#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DEBUG=0

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))

cd pythia8

pidArr=()
NAMES=""
for (( i=1; i<=$NUM_PROC; i++ ))
do
    P8FILE=$(python settings.py $NUM_EVT $JOB_TYPE $NUM_PROC $i)
    ./pythia8.exe $JOB_TYPE $P8FILE &
    pidArr+=($!)
    pidArr+=" "
    NAMES+="particles_pythia8_"$P8FILE".root"
    NAMES+=" "
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done

MERGE="particles_pythia8_"$(python -c "import sys; word = sys.argv[1]; print word[0:-2]" $P8FILE)".root"

if [ $NUM_PROC -gt 1 ]; then
    hadd -f $MERGE $NAMES

    for tmp in $NAMES
    do
        rm $tmp
    done
else
    mv $NAMES $MERGE
fi

if [ $DEBUG -eq 0 ]; then
    REMAIN=$(python -c "import sys; word = sys.argv[1]; print word[0:-2]" $P8FILE)
    rm $REMAIN*.cmnd
fi

mv $MERGE ..

exit
