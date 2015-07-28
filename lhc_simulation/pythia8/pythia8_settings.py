#!/usr/bin/python

import sys
import os

# Hard coded tune type
# 0: 4C from Herwig++ defaults (CTEQ6L1)
# 1: CUETP8S1 from CMS         (CTEQ6L1)
tune = 0

# Hard coded PDF choice
# 0: CTEQ6L1
# 1: CT10
# 2: MSTW2008LO
pdf = 0

# Read run settings from command line parameters
if len(sys.argv) != 5:
    print "Incorrect amount of command line parameters"
    sys.exit()

tot_evts = int(sys.argv[1])
mode = int(sys.argv[2])
procs = int(sys.argv[3])
proc_id = int(sys.argv[4])

if tot_evts%procs!=0:
    print "Number of processors and events should fit"
    sys.exit()

# Set a unique name for the run
name = ""
if mode == 1:
    name += "dijet"
elif mode == 2:
    name += "gammajet"
elif mode == 3:
    name += "Zjet"
elif mode == 4:
    name += "ttbar"

name += "_"
name += str(tot_evts)

name += "_"
name += str(proc_id)

print name
name += ".cmnd"

f = open(name,'w')

f.write("Random:setSeed = on\n")
f.write("Random:seed = {}\n\n".format(proc_id*10000) )

f.write("Main:numberOfEvents = {}\n\n".format(tot_evts/procs))

f.write("! Generic settings\n")
f.write("Next:numberShowInfo = 0\n")
f.write("Next:numberShowProcess = 0\n")
f.write("Next:numberShowEvent = 0\n")
f.write("Next:numberCount = 0\n\n")

f.write("! Set particles with long enough lifetimes to stable and photon radiation\n")
f.write("ParticleDecays:allowPhotonRadiation = on\n");
f.write("ParticleDecays:limitTau0=on\n")
f.write("ParticleDecays:tau0Max=10.\n\n")

f.write("! Event weighting\n")
f.write("PhaseSpace:bias2Selection = on\n")
f.write("PhaseSpace:bias2SelectionPow = 4.5\n")
f.write("PhaseSpace:bias2SelectionRef = 15.\n\n")

f.write("! CM energy\n")
f.write("Beams:eCM = 8000.\n\n")

if mode==1:
    f.write("HardQCD:all = on\n")
    f.write("PhaseSpace:pTHatMin = 30.\n\n")
if mode==2:
    f.write("PromptPhoton:qg2qgamma = on\n");
    f.write("PromptPhoton:qqbar2ggamma = on\n");
    f.write("PromptPhoton:gg2ggamma = on\n");
    f.write("PhaseSpace:pTHatMin = 10.\n\n");
if mode==3:
    f.write("! Produce both Z0's and gammas\n");
    f.write("WeakZ0:gmZmode = 0\n");
    f.write("WeakBosonAndParton:qqbar2gmZg = on\n");
    f.write("WeakBosonAndParton:qg2gmZq  = on\n");
    f.write("23:onMode = off\n");
    f.write("23:7:onMode = on\n");
    f.write("PhaseSpace:pTHatMin = 20.\n\n");
    f.write("PhaseSpace:mHatMin = 75.\n");
if mode==4:
    f.write("Top::gg2ttbar = on\n");
    f.write("Top::qqbar2ttbar = on\n");
    f.write("PhaseSpace:pTHatMin = 30.\n\n")

f.write("Tune:preferLHAPDF = 2\n")
if tune==0:
    f.write("! Tune (4C)\n")
    f.write("Tune:ee = 3\n")
    f.write("Tune:pp = 5\n")
elif tune==1:
    f.write("! CMS UE Tune CUETP8S1-CTEQ6L1\n")
    #f.write("Tune:pp = 15\n\n")
    f.write('Tune:ee 3\n')
    f.write('Tune:pp 5\n')
    f.write('MultipartonInteractions:pT0Ref=2.1006\n')
    f.write('MultipartonInteractions:ecmPow=0.21057\n')
    f.write('MultipartonInteractions:expPow=1.6089\n')
    f.write('MultipartonInteractions:a1=0.00\n')
    f.write('ColourReconnection:range=3.31257\n\n')

if pdf==0:
    f.write("PDF:pSet = LHAPDF6:cteq6l1\n\n")

#f.write("PartonLevel:MPI = off\n")
#f.write("PartonLevel:ISR = off\n")
#f.write("PartonLevel:FSR = off\n")

#f.write("HadronLevel:Hadronize = off\n")

f.close()
