#!/bin/bash

source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.00.00/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

export LD_LIBRARY_PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/lib:/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/lib:/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/lib64:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/bin:/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/bin:$PATH

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DEBUG=0

WRKDIR=/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/herwig

pidArr=()
NAMES=""
for (( i=1; i<=$NUM_PROC; i++ ))
do
    cp $WRKDIR/settings.py .
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

xrdcp $MERGE root://eoscms.cern.ch//eos/cms/store/group/phys_jetmet/hsiikone/.

exit
