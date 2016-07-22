#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
OFFSET_IDX=$4
DEBUG=0
BODY=""

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))

cd pythia8

pidArr=()
NAMES=""
for (( i=0; i<$NUM_PROC; i++ ))
do
    P8FILE=$(python settings.py $NUM_EVT $JOB_TYPE $NUM_PROC $(($i+$OFFSET_IDX)))
    echo $JOB_TYPE
    echo $P8FILE
    ./pythia8.exe $JOB_TYPE $P8FILE &
    pidArr+=($!)
    pidArr+=" "
    NAMES+="particles_pythia8_"$P8FILE".root"
    NAMES+=" "
    if [ $i -eq 0 ]; then
        BODY+=$(python -c "import sys; import re; word = sys.argv[1]; print re.search('^(.+?)_[0-9]+$',word).group(1)" $P8FILE)
    fi
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done
MERGE="particles_pythia8_"$BODY".root"
echo $NAMES
echo $MERGE
echo $NUM_PROC

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
    rm $BODY*.cmnd
fi

mv $MERGE ..

exit
