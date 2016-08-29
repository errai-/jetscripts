#ifndef JETSORTER_H
#define JETSORTER_H

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

#include <cmath>
#include <ctime>
// FastJet interface
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/PseudoJet.hh"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"
// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

// tdrStyle
#include "../src/tdrstyle_mod1.C"
// scripts
#include "../../JetSorter/jetsorter_auxiliary.h"
#include "../include/help_functions.h"
#include "../include/MinimalEvent.h"

class JetSorter{
private:  

  // ROOT initialization:  
  static const int ptBins = 48.;
  const double ptRange[ptBins+1]=
  //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
  //3637, 3832, 4037};//  
  
  // Create the ROOT application environment.
  //TApplication* theApp;
  TFile* outFile;
  TH1D* ptProfile;
  TH1D* jetMultipl;
  TProfile* gluonQuark;
  vector<TH1D*> chargeIndicator;
  
  // Book histograms.
  vector<TProfile*> fractionProfilesGluon;
  vector<TProfile*> fractionProfilesQuark;
  vector<TProfile*> fractionProfilesLQuark;
  vector<TProfile*> fractionProfilesHQuark;
  vector<TProfile*> fractionProfilesAll;
  // ROOT initialization ^
  
  // Fastjet initialization:
  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range

  // Fastjet analysis - select algorithm and parameters
  fastjet::JetDefinition* jetDef;

  // Fastjet initialization ^
  
  // Event loop initialization:
  MinimalEvent eventHandler;
  std::ifstream input;
  size_t iEvent = 0;

  // Initialize energy counters 
  double piPlus, piMinus,  pi0Gamma, gamma, 
    kaPlus, kaMinus, kSZero, kLZero, 
    proton, aproton, neutron, aneutron,
    lambda0, sigma, elecmuon,
    others, etSum;

  vector<fastjet::PseudoJet> sortedJets; 
  vector<fastjet::PseudoJet> jetParts;

  vector<fastjet::PseudoJet> fjInputs;
public:
  
  JetSorter(){
    outFile = new TFile("sortedjets.root", "RECREATE");
    ptProfile = new TH1D("pT bins","Pt bins", ptBins, ptRange);
    jetMultipl = new TH1D("Jet multiplicity","Jet multiplicity",50,0.,50.);
    gluonQuark = new TProfile("gq","gq",ptBins,ptRange);
    
    InitCI();
    
    InitFP();
    
    jetDef = new fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
    
    input.open("pythia8data.txt", std::ifstream::in);
  }

  ~JetSorter(){
    delete outFile;
    delete ptProfile;
    delete jetMultipl;
    delete gluonQuark;
    delete jetDef;
    
    for (auto it = chargeIndicator.begin(); it != chargeIndicator.end(); ++it){ delete (*it);}
    for (auto it = fractionProfilesAll.begin(); it != fractionProfilesAll.end(); ++it){ delete (*it);}
    for (auto it = fractionProfilesHQuark.begin(); it != fractionProfilesHQuark.end(); ++it){ delete (*it);}
    for (auto it = fractionProfilesLQuark.begin(); it != fractionProfilesLQuark.end(); ++it){ delete (*it);}
    for (auto it = fractionProfilesGluon.begin(); it != fractionProfilesGluon.end(); ++it){ delete (*it);}
    for (auto it = fractionProfilesQuark.begin(); it != fractionProfilesQuark.end(); ++it){ delete (*it);}
  }

  void InitCI();

  void InitFP();

  void EventLoop();
  
  bool ParticlesToJetsorterInput();
  
  void JetLoop();
  
  void ParticleLoop();
  
  void HistFill(int);

  void WriteResults();

};



#endif // JETSORTER_H