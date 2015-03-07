
#ifndef JETANALYSIS_H
#define JETANALYSIS_H

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
#include "SimEvent.h"
#include "JetEvent.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include "help_functions.h"
// tdrStyle
#include "tdrstyle_mod1.C"

using std::cout;
using std::endl;

class JetAnalysis {
private: 
// Mostly Fastjet and histogramming
  
///////////////////////
// ROOT initialization:
///////////////////////
  static const int ptBins = 48.;
  const double ptRange[ptBins+1]=
  //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
  //3637, 3832, 4037};//  
  
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
////////////////////////
// ROOT initialization ^
////////////////////////
  
//////////////////////////
// Fastjet initialization:
//////////////////////////
  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range

  // Analysis: select algorithm and parameters
  fastjet::JetDefinition* jetDef;

  // PseudoJet storages
  vector<fastjet::PseudoJet> sortedJets, sortedGhosts, jetParts, fjInputs; 
///////////////////////////
// Fastjet initialization ^
///////////////////////////
    
  Timer timer;
  
  // Event loop:
  std::ifstream input;
  size_t iEvent = 0;

  // Energy counters: 
  double piPlus, piMinus,  pi0Gamma, gamma, 
    kaPlus, kaMinus, kSZero, kLZero, 
    proton, aproton, neutron, aneutron,
    lambda0, sigma, elecmuon,
    others, etSum;
  
  // Jet flavor studies: 
  vector<size_t> FlavorIndices; 
  int partonHadronFlavour[2];
  
  // Weights etc.:
  double partSum, chargSum, chargWSum, chargW2Sum, w2;
  int quarkJetCharge, isHadron;
  
public : 
// Interface to the tree and functions

/////////
// Input:
/////////  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  // If this is too small, Segfaults may follow.
  static const Int_t kMaxfParts = 5000;
  // Declaration of leaf types
  //SimEvent        *event;
  Int_t           fParts_;
  Double_t        fParts_fPx[kMaxfParts];   //[fParts_]
  Double_t        fParts_fPy[kMaxfParts];   //[fParts_]
  Double_t        fParts_fPz[kMaxfParts];   //[fParts_]
  Double_t        fParts_fE[kMaxfParts];   //[fParts_]
  Int_t           fParts_fPDGCode[kMaxfParts];   //[fParts_]
  Int_t           fParts_fChargeTimes3[kMaxfParts];   //[fParts_]
  Bool_t          fParts_IsPi0Photon[kMaxfParts];   //[fParts_]
  Bool_t          fParts_IsJetFlavor[kMaxfParts];   //[fParts_]
  Bool_t          fParts_IsExcitedState[kMaxfParts];   //[fParts_]

  // List of branches
  TBranch        *b_event_fParts_;   //!
  TBranch        *b_fParts_fPx;   //!
  TBranch        *b_fParts_fPy;   //!
  TBranch        *b_fParts_fPz;   //!
  TBranch        *b_fParts_fE;   //!
  TBranch        *b_fParts_fPDGCode;   //!
  TBranch        *b_fParts_fChargeTimes3;   //!
  TBranch        *b_fParts_IsPi0Photon;   //!
  TBranch        *b_fParts_IsJetFlavor;   //!
  TBranch        *b_fParts_IsExcitedState;   //!
//////////
// Input ^
//////////

  // Output:
  TTree          *outTree;
  JetEvent jEvent;
  TBranch jetBranch;
  
  JetAnalysis(TTree *tree = 0) : fChain(0) {
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
  
    outTree = new TTree("JetTree","Tree filled with simulated jet data.");
   
    // Create a ROOT Tree and one superbranch
    outTree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
    outTree->SetCacheSize(10000000);  // set a 10 MBytes cache (useless when writing local files)
    TTree::SetBranchStyle(1);
  }

  ~JetAnalysis() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
    
    delete outFile;
    delete jetDef;
  }
   
  virtual Int_t    GetEntry(Long64_t);
  virtual Long64_t LoadTree(Long64_t);
  virtual void     Init(TTree *tree); // Chain
  virtual void     InitCI(); // Charge indicator
  virtual void     InitFP(); // Fraction profiles
  virtual void     Show(Long64_t = -1);
//   virtual void     Insert(double,double=0,double=0,double=0,int=0,int=0);
//   virtual void     Clear();
  virtual void     EventLoop();
  virtual bool     ParticlesToJetsorterInput();
  virtual void     JetLoop();
  virtual void     FlavorLoop(size_t);
  virtual void     ParticleLoop(size_t);
  virtual void     HistFill(int);
  virtual void     FillerHandle(vector<TProfile*> &, double, double, double,
    double, double, double, double, double, double, double, double, double, 
    double, double, double, double, double, double);
  virtual void     WriteResults();
  
};

#endif // JETANALYSIS_H