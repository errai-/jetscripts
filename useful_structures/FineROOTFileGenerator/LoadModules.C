#include <iostream>
#include <TROOT.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <vector>

//#include <QCDModules/QCDPFJet.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void LoadModules()
{
  gROOT->ProcessLine(".L QCDModules/QCDEvent.cc++");
  gROOT->ProcessLine(".L QCDModules/QCDJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDMET.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDPFJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDEventHdr.cc++");
  gROOT->ProcessLine(".L QCDModules/LorentzVector.h++");
  gROOT->ProcessLine(".L Auxiliary.C++");
}

