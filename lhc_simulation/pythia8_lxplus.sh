#!/bin/bash


export PYTHIA8=/afs/cern.ch/user/h/hsiikone/Cern/installs/include                                       
export PYTHIA8DATA=/afs/cern.ch/user/h/hsiikone/Cern/installs/share/Pythia8/xmldoc                      
                                                                                                        
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.00/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh         
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh                                         
                                                                                                        
export LD_LIBRARY_PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/lib:/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/lib:/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/lib64:$LD_LIBRARY_PATH                                                                                    
export PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/bin:/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/bin:$PATH                                                                            

NUM_EVT=$1
JOB_TYPE=$2
NUM_PROC=$3
DEBUG=0

EVT_PER_RUN=$(($NUM_EVT/$NUM_PROC))
WRKDIR=/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/pythia8

pidArr=()
NAMES=""
for (( i=1; i<=$NUM_PROC; i++ ))
do
    P8FILE=$(python $WRKDIR/pythia8_settings.py $NUM_EVT $JOB_TYPE $NUM_PROC $i)
    $WRKDIR/pythia8.exe $JOB_TYPE $P8FILE &
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

xrdcp $MERGE root://eoscms.cern.ch//eos/cms/store/group/phys_jetmet/hsiikone/.
rm $MERGE

exit

