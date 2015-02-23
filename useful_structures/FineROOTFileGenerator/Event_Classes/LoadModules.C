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

//#include <QCDModules/QCDPFJet.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void LoadModules()
{
  gROOT->ProcessLine(".L QCDModules/QCDJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDMET.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDPFJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDEventHdr.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDEvent.cc+");
  //gROOT->ProcessLine(".L Auxiliary.C++");
}

