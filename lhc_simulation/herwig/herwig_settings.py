#!/usr/bin/python

import sys
import os

# Hard coded tune type
# 0: EE 5C from Herwig++ defaults (CTEQ6L1)
# 1: CUETHS1 from CMS             (CTEQ6L1)
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

name += "_"
name += str(tot_evts)

name += "_"
name += str(proc_id)

name += ".in"
print name

f = open(name,'w')

f.write('################################################################\n')
f.write('# This file is based on the CUETHS1 Herwig++ tune from CMSSW,\n')
f.write('# which is based on the generic Herwig++ tune EE-5C for CTEQ6l1.\n')
f.write('# See CMSSW/Configuration/Generator/python Hpp configurations.\n')
f.write('################################################################\n\n')

f.write('################################\n')
f.write('# Setting up the event generator\n')
f.write('################################\n')
f.write('cd /Herwig/Generators\n')
f.write('set LHCGenerator:RandomNumberGenerator:Seed {}\n'.format(100*proc_id))
f.write('set LHCGenerator:NumberOfEvents {}\n'.format(int(tot_evts/procs)))
f.write('set LHCGenerator:DebugLevel 1\n')
f.write('set LHCGenerator:UseStdout 0\n')
f.write('set LHCGenerator:PrintEvent 0\n')
f.write('set LHCGenerator:MaxErrors 10000\n\n')

f.write('# Make some particles stable, according to their nominal lifetimes\n')
f.write('set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n\n')

f.write('########################\n')
f.write('# LHC running parameters\n')
f.write('########################\n\n')

f.write('# Energy extrapolation\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV\n')
if tune==0:
    f.write('set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33\n')
elif tune==1:
    f.write('set /Herwig/UnderlyingEvent/MPIHandler:Power 0.3705288\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV\n\n')

f.write('# Colour reconnection settings\n')
f.write('set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes\n')
if tune==0:
    f.write('set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49\n\n')
elif tune==1:
    f.write('set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.5278926\n\n')

f.write('# Colour Disrupt settings\n')
if tune==0:
    f.write('set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80\n\n')
elif tune==1:
    f.write('set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.6284222\n\n')

f.write('# inverse hadron radius\n')
if tune==0:
    f.write('set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30\n\n')
elif tune==1:
    f.write('set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.254998\n\n')

f.write('# MPI model settings\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes\n')
f.write('set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2\n\n')

f.write('# LHAPDF settings\n')
f.write('cd /Herwig/Partons\n')
f.write('create ThePEG::LHAPDF customPDF ThePEGLHAPDF.so\n')
if pdf==0:
    f.write('set customPDF:PDFName cteq6l1\n')
elif pdf==1:
    f.write('set customPDF:PDFName CT10\n')
f.write('set customPDF:RemnantHandler /Herwig/Partons/HadronRemnants\n')
f.write('set /Herwig/Particles/p+:PDF customPDF\n')
f.write('set /Herwig/Particles/pbar-:PDF customPDF\n\n')

f.write('# CM energy\n')
f.write('set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 8000.0*GeV\n\n')

f.write('# Intrinsic pT tune extrapolated to LHC energy\n')
f.write('set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.0*GeV\n\n')

f.write('# ptHat min\n')
if mode==1:
    f.write('set /Herwig/Cuts/JetKtCut:MinKT 30.0*GeV\n\n')

f.write('##############################################\n')
f.write('# Matrix Elements for hadron-hadron collisions\n')
f.write('##############################################\n\n')

f.write('# Event weighting scheme\n')
f.write('mkdir /Herwig/Weights\n')
f.write('cd /Herwig/Weights\n')
f.write('create ThePEG::ReweightMinPT reweightMinPT ReweightMinPT.so\n')
f.write('set reweightMinPT:Power 4.5\n')
f.write('set reweightMinPT:Scale 15*GeV\n\n')

f.write('# Set matrix element settings\n')
f.write('cd /Herwig/MatrixElements/\n')
if mode==1:
    f.write('insert SimpleQCD:MatrixElements[0] MEQCD2to2\n')
f.write('insert SimpleQCD:Preweights[0] /Herwig/Weights/reweightMinPT\n\n')

f.write('# Save final particles and hardest subprocess particles\n')
f.write('cd /Herwig/Analysis\n')
f.write('create jetanalysis::StoreParticles jetAnalysis ../lib/libStoreParticles.so\n')
f.write('insert /Herwig/Generators/LHCGenerator:AnalysisHandlers 0 jetAnalysis\n\n')

f.write('# For now saverun does not work with LHAPDF\n')
f.write('run {} /Herwig/Generators/LHCGenerator\n'.format(name[0:-3]))

f.close()