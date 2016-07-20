#!/bin/bash

export PYTHIA8=/afs/cern.ch/user/h/hsiikone/Cern/installs/include
export PYTHIA8DATA=/afs/cern.ch/user/h/hsiikone/Cern/installs/share/Pythia8/xmldoc

source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.00.00/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

export LD_LIBRARY_PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/lib:/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/lib:/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/lib64:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/bin:/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/bin:$PATH

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DIVINE_IDX=$4
DEBUG=0
BODY=""

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))
WRKDIR=/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/pythia8

pidArr=()
NAMES=""
for (( i=0; i<$NUM_PROC; i++ ))
do
    cp $WRKDIR/settings.py .
    P8FILE=$(python settings.py $NUM_EVT $JOB_TYPE $NUM_PROC $(($i+$DIVINE_IDX)))
    $WRKDIR/pythia8.exe $JOB_TYPE $P8FILE &
    pidArr+=($!)
    pidArr+=" "
    NAMES+="particles_pythia8_"$P8FILE".root"
    NAMES+=" "
    if [ $i -eq 0 ]; then
        BODY+=$P8FILE
    fi
done

for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
done

MERGE="particles_pythia8_"$BODY"f.root"

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
    rm $(python -c "import sys; import re; word = sys.argv[1]; print re.search('^(.+?)_[0-9]+$',word).group(1)" $BODY)*.cmnd
fi

xrdcp $MERGE root://eoscms.cern.ch//eos/cms/store/group/phys_jetmet/hsiikone/.
rm $MERGE

exit

