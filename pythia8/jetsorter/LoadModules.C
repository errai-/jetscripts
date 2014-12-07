#include <iostream>
#include <TROOT.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <vector>

//#include <Hannouris/QCDPFJet.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void LoadModules()
{
  gROOT->ProcessLine(".L Hannouris/QCDEvent.cc++");
  gROOT->ProcessLine(".L Hannouris/QCDJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDMET.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDPFJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDEventHdr.cc++");
  gROOT->ProcessLine(".L Hannouris/LorentzVector.h++");
  gROOT->ProcessLine(".L Auxiliary.C++");
}

