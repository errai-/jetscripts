#!/bin/bash

setenv RFIO_USE_CASTOR_V2 yes

setenv STAGE_HOST castorpublic

setenv STAGE_SVCCLASS default

source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

export LD_LIBRARY_PATH=/afs/cern.ch/user/h/hsiikone/Cern/pythia8/lib:/afs/cern.ch/user/h/hsiikone/Cern/fastjet/lib:/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.9.2/x86_64-slc6-gcc49-opt/lib64:$LD_LIBRARY_PATH

export PATH=/afs/cern.ch/user/h/hsiikone/Cern/pythia8/bin:/afs/cern.ch/user/h/hsiikone/Cern/fastjet/bin:/afs/cern.ch/sw/lcg/contrib/gcc/4.9.2/x86_64-slc6-gcc49-opt/bin:$PATH

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=1

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))

pidArr=()
for (( i=1; i<=$NUM_PROC; i++ ))
do
    /afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/pythia8/pythia8.exe $EVT_PER_RUN $JOB_TYPE $i $NUM_PROC &
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

xrdcp $MERGE root://eoscms.cern.ch//eos/cms/store/group/phys_jetmet/hsiikone/.

exit
