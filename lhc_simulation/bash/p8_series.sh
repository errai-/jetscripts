#!/bin/bash

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=8
DIVINE_IDX=0
NUM_LAYERS=4

EVT_PER_LAYER=$(($NUM_EVT/$NUM_LAYERS))

for (( i=1; i<=$NUM_LAYERS; i++ ))
do
    bsub -n $NUM_PROC -q 1nw -R "pool>20000" pythia8_lxplus.sh $EVT_PER_LAYER $JOB_TYPE $NUM_PROC
done
