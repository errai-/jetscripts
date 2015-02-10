
#ifndef ROOTJETSORT_H
#define ROOTJETSORT_H

// General
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <THStack.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

// FastJet interface
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// Header file for the classes stored in the TTree if any.
#include "MinimalEvent.h"
// Fixed size dimensions of array or collections stored in the TTree if any.
#include "help_functions.h"
// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
//#include "../../JetSorter/jetsorter_auxiliary.h"

using std::cout;
using std::endl;

class RootJetSort {
private: // Mostly Fastjet and histogramming
  
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
  
public : // Interface to the tree and functions
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  //MinimalEvent    *event;
  vector<double>  px;
  vector<double>  py;
  vector<double>  pz;
  vector<double>  e;
  vector<int>     id;
  ULong_t         particles;
  Timer timer;

  // List of branches
  TBranch        *b_event_px;   //!
  TBranch        *b_event_py;   //!
  TBranch        *b_event_pz;   //!
  TBranch        *b_event_e;   //!
  TBranch        *b_event_id;   //!
  TBranch        *b_event_particles;   //!

  RootJetSort(TTree *tree = 0) : fChain(0) {
    if (tree == 0) {
       TChain * chain = new TChain("Pythia8Tree","");
      chain->Add("particle_storage.root/Pythia8Tree;1");
      tree = chain;
    }
    Init(tree);
    
    outFile = new TFile("sortedjets.root", "RECREATE");
    ptProfile = new TH1D("pT bins","Pt bins", ptBins, ptRange);
    jetMultipl = new TH1D("Jet multiplicity","Jet multiplicity",50,0.,50.);
    gluonQuark = new TProfile("gq","gq",ptBins,ptRange);
    
    InitCI();
    InitFP();
    
    jetDef = new fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
  }

  ~RootJetSort() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
    
    delete outFile;
//     delete ptProfile;
//     delete jetMultipl;
//     delete gluonQuark;
    delete jetDef;
    
//     for (auto it = chargeIndicator.begin(); it != chargeIndicator.end(); ++it){ delete (*it);}
//     for (auto it = fractionProfilesAll.begin(); it != fractionProfilesAll.end(); ++it){ delete (*it);}
//     for (auto it = fractionProfilesHQuark.begin(); it != fractionProfilesHQuark.end(); ++it){ delete (*it);}
//     for (auto it = fractionProfilesLQuark.begin(); it != fractionProfilesLQuark.end(); ++it){ delete (*it);}
//     for (auto it = fractionProfilesGluon.begin(); it != fractionProfilesGluon.end(); ++it){ delete (*it);}
//     for (auto it = fractionProfilesQuark.begin(); it != fractionProfilesQuark.end(); ++it){ delete (*it);}
 }
   
  virtual Int_t    GetEntry(Long64_t);
  virtual Long64_t LoadTree(Long64_t);
  virtual void     Init(TTree *tree); // Chain
  virtual void     InitCI(); // Charge indicator
  virtual void     InitFP(); // Fraction profiles
  virtual void     Show(Long64_t = -1);
  virtual void     Insert(double,double=0,double=0,double=0,int=0,int=0);
  virtual void     Clear();
  virtual void     EventLoop();
  virtual bool     ParticlesToJetsorterInput();
  virtual void     JetLoop();
  virtual void     ParticleLoop();
  virtual void     HistFill(int);
  virtual void     WriteResults();
  
  vector<fastjet::PseudoJet> Momentum;
  vector<int> status;
  vector<int> id_store;   
};

#endif // ROOTJETSORT_H


