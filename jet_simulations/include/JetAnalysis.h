/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 7.3.2015                     //
/////////////////////////////////////////////////

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
#include <TLorentzVector.h>

// FastJet interface
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// Header file for the classes stored in the TTree if any.
#include "PrtclEvent.h"
#include "JetEvent.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include "help_functions.h"
// tdrStyle
#include "tdrstyle_mod1.C"

using std::cout;
using std::endl;
using std::string;

class JetAnalysis {
private: 
  
//////////////////////////
// Fastjet initialization:
//////////////////////////
  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range
  int jetsPerEvent = 2;   // How many leading jets are stored in a run

  // Analysis: select algorithm and parameters
  fastjet::JetDefinition* jetDef;

  // PseudoJet storages
  vector<fastjet::PseudoJet> sortedJets, sortedGhosts, jetParts, fjInputs; 

//////////
// Others:
//////////
  Timer timer;
  
  // Event loop:
  std::ifstream input;
  size_t iEvent = 0;

  // Energy counters: 
  TLorentzVector piPlus, piMinus,  pi0Gamma, gamma, 
    kaPlus, kaMinus, kSZero, kLZero, 
    proton, aproton, neutron, aneutron,
    lambda0, sigma, elec, muon,
    others, etSum;
  
  // Jet flavor studies: 
  vector<size_t> FlavorIndices; 
  
  // Weights etc.:
  double partSum, chargSum, chargWSum, chargW2Sum, w2;
  int quarkJetCharge, flavour;
  
  // Temporary fractions:
  double chf;
  double nhf;
  double phf;
  double elf;
  double muf;
  
  // Temporary masses:
  double chm;
  double nhm;
  double phm;
  double elm;
  double mum;
  
public : 

/////////
// Input:
/////////  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  // If this is too small, Segfaults may follow.
  static const Int_t kMaxfParts = 5000;
  // Declaration of leaf types
  //PrtclEvent        *event;
  Int_t           fParts_;
  Double_t        fParts_fP4_fCoordinates_fX[kMaxfParts];   //[fParts_]
  Double_t        fParts_fP4_fCoordinates_fY[kMaxfParts];   //[fParts_]
  Double_t        fParts_fP4_fCoordinates_fZ[kMaxfParts];   //[fParts_]
  Double_t        fParts_fP4_fCoordinates_fT[kMaxfParts];   //[fParts_]
  
  Int_t           fParts_fPDGCode[kMaxfParts];   //[fParts_]
  Int_t           fParts_fChargeTimes3[kMaxfParts];   //[fParts_]
  
  Bool_t          fParts_IsPi0Photon[kMaxfParts];   //[fParts_]
  Bool_t          fParts_IsJetFlavor[kMaxfParts];   //[fParts_]
  Bool_t          fParts_IsExcitedState[kMaxfParts];   //[fParts_]

  // List of branches
  TBranch        *b_event_fParts_;   //!
  TBranch        *b_fParts_fP4_fCoordinates_fX;   //!
  TBranch        *b_fParts_fP4_fCoordinates_fY;   //!
  TBranch        *b_fParts_fP4_fCoordinates_fZ;   //!
  TBranch        *b_fParts_fP4_fCoordinates_fT;   //!
  
  TBranch        *b_fParts_fPDGCode;   //!
  TBranch        *b_fParts_fChargeTimes3;   //!
  
  TBranch        *b_fParts_IsPi0Photon;   //!
  TBranch        *b_fParts_IsJetFlavor;   //!
  TBranch        *b_fParts_IsExcitedState;   //!
//////////
// Output:
//////////
  TFile          *outFile;
  TTree          *outTree;
  TBranch        *jetBranch;
  JetEvent       *jEvent;

/////////////
// Functions:
/////////////
  JetAnalysis(TTree * = 0);

  ~JetAnalysis();
   
  virtual void     Init(TTree *tree); // Chain
  
  virtual Int_t    GetEntry(Long64_t);
  virtual Long64_t LoadTree(Long64_t);
  virtual void     Show(Long64_t = -1);
  
  virtual void     EventLoop();
  virtual bool     ParticlesToJetsorterInput();
  virtual void     JetLoop();
  virtual void     FlavorLoop(size_t);
  virtual void     ParticleLoop(size_t);
  virtual void     TypeSort();
};

#endif // JETANALYSIS_H