#!/bin/bash

source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.00.00/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

export LD_LIBRARY_PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/lib:/afs/cern.ch/user/h/hsiikone/Cern/jetscripts/lhc_simulation/lib:/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/lib64:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/user/h/hsiikone/Cern/installs/bin:/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/bin:$PATH

FILE=$1
FINALTAG=$2
SUBDIR=$3

WRKDIR=store/group/phys_jetmet/hsiikone/$SUBDIR
TOTALPATH=root://eoscms.cern.ch//eos/cms/$WRKDIR

COMMAND=$TOTALPATH/$FILE
COMMAND+="_"
COMMAND+=$FINALTAG
COMMAND+=".root"
COMBINE=$(ls $TOTALPATH/$FILE*)
for fname in `cat $COMBINE`; do
    $COMMAND+=" "
    $COMMAND+=$TOTALPATH/$fname
done
echo $COMMAND

hadd -f $COMMAND

exit
