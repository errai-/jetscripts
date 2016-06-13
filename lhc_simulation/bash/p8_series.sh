#!/bin/bash

NUM_EVT=$1
NUM_LAYERS=$2
JOB_TYPE=$3
NUM_PROC=8

EVT_PER_LAYER=$(($NUM_EVT/$NUM_LAYERS))

for (( i=0; i<$NUM_LAYERS; i++ ))
do
    bsub -n $NUM_PROC -q 2nw -R "pool>20000" pythia8_lxplus.sh $EVT_PER_LAYER $JOB_TYPE $NUM_PROC $(($i*$NUM_PROC))
done
