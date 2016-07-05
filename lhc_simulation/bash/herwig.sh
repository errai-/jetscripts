#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DEBUG=1

cd herwig

pidArr=()
NAMES=""
for (( i=1; i<=$NUM_PROC; i++ ))
do
    HFILE=$(python settings.py $NUM_EVT $JOB_TYPE $NUM_PROC $i)
    Herwig read $HFILE &
    pidArr+=($!)
    pidArr+=" "
    NAMES+="particles_herwig_"$(python -c "import sys; word = sys.argv[1]; print word[0:-3]" $HFILE)".root"
    NAMES+=" "
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done

MERGE="particles_herwig_"$(python -c "import sys; word = sys.argv[1]; print word[0:-5]" $HFILE)".root"

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
    REMAIN=$(python -c "import sys; word = sys.argv[1]; print word[0:-5]" $HFILE)
    rm $REMAIN*.in
    rm $REMAIN*.out
    rm $REMAIN*.log
    rm $REMAIN*.tex
fi

mv $MERGE ..

exit
