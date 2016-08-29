#include <iostream>
#include <TROOT.h>
#include <TMath.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzE4D.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void LoadModules()
{
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDJet.cc+");
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDMET.cc+");
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDPFJet.cc+");
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDEventHdr.cc+");
  gROOT->ProcessLine(".L KKousour/QCDAnalysis/src/QCDEvent.cc+");
}

