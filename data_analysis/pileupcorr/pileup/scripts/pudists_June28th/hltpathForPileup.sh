#!/bin/bash
##############################################################################

# This script will produce customized pile-up profiles for 2012 jet triggers #

# Update lumiSummary*.json and pileup_JSON_*.txt below as needed             #

##############################################################################
JSON=lumiSummary_June28th_3fb.json
PUJSON=pileup_JSON_DCSONLY_190389-196531_patch2.txt



# Get the centrally produced pile-up JSON file:

# https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#2012_Pileup_JSON_Files

# rfdir /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/

# rfcp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/$PUJSON .



# Customize pile-up JSON for each trigger

# Replace lumiCalc2.py below with pixelLumiCalc.py, when available for 2012

# https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Modifying_Pileup_JSON_for_Use_wi



#echo "Using input JSON file: "$JSON

#echo "Using pile-up JSON file: "$PUJSON



echo "HLT_PFJet40_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet40_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet40_*" -o lumicalc2_jt40.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt40.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt40.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt40.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt40 pileup12_jt40.root --verbose



echo "HLT_PFJet80_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet80_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet80_*" -o lumicalc2_jt80.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt80.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt80.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt80.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt80 pileup12_jt80.root --verbose



echo "HLT_PFJet140_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet140_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet140_*" -o lumicalc2_jt140.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt140.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt140.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt140.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt140 pileup12_jt140.root --verbose



echo "HLT_PFJet200_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet200_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet200_*" -o lumicalc2_jt200.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt200.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt200.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt200.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt200 pileup12_jt200.root --verbose



echo "HLT_PFJet260_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet260_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet260_*" -o lumicalc2_jt260.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt260.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt260.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt260.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt260 pileup12_jt260.root --verbose



echo "HLT_PFJet320_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet320_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet320_*" -o lumicalc2_jt320.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt320.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt320.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt320.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt320 pileup12_jt320.root --verbose



echo "HLT_PFJet400_*"

#lumiCalc2.py recorded -i $JSON --hltpath "HLT_PFJet400_*"

lumiCalc2.py lumibyls -i $JSON --hltpath "HLT_PFJet400_*" -o lumicalc2_jt400.csv

pileupReCalc_HLTpaths.py -i lumicalc2_jt400.csv --inputLumiJSON $PUJSON -o pileupJSON12_jt400.txt

pileupCalc.py -i $JSON --inputLumiJSON=pileupJSON12_jt400.txt --calcMode=true --minBiasXsec=71000 --maxPileupBin=60 --numPileupBins=600 --pileupHistName=pileup_jt400 pileup12_jt400.root --verbose







echo "Adding pileup12_*.root files together in pileup12_May16th.root"

hadd pileup12_June28th_mb710.root pileup12_jt40.root pileup12_jt80.root pileup12_jt140.root pileup12_jt200.root pileup12_jt260.root pileup12_jt320.root pileup12_jt400.root 
