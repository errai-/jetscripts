/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 12.3.2015                    //
/////////////////////////////////////////////////

#ifndef JETANALYSIS_H
#define JETANALYSIS_H

// General
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

// FastJet interface
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// Header file for the classes stored in the TTree if any.
#include "events/PrtclEvent.h"
#include "events/JetEvent.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include "generic/help_functions.h"

using std::cout;
using std::endl;
using std::vector;

class JetAnalysis 
{
public : 
   JetAnalysis(TTree * = 0);

   ~JetAnalysis();
   
   virtual void     Init(TTree*, char* = "jet_storage.root"); // Chain
  
   virtual Int_t    GetEntry(Long64_t);
   virtual Long64_t LoadTree(Long64_t);
   virtual void     Show(Long64_t = -1);
  
   virtual void     EventLoop();
   virtual void     ParticlesToJetsorterInput();
   virtual void     GhostsToJetsorterInput();
   virtual void     JetLoop();
   virtual void     FlavorLoop(size_t);
   virtual void     ParticleLoop(size_t);
   virtual void     TypeSort();

private: 
///////////
// Fastjet:
///////////
   int power     = -1;     /* -1 = ant-kT; 0 = C/A; 1 = kT */
   double R      = 0.5;    /* Jet size. */
   double pTMin  = 20.0;   /* Min jet pT */
   double etaMax = 1.3;    /* Pseudorapidity range */
   int jetsPerEvent = 2;   /* How many leading jets are stored in a run /

   /* Analysis: select algorithm and parameters */
   fastjet::JetDefinition* jetDef;

   /* PseudoJet storages */
   vector<fastjet::PseudoJet> sortedJets, sortedGhosts, jetParts, fjInputs; 

/////////
// Input:
/////////  
   TTree          *fChain;   //! 
   Int_t           fCurrent; //! current Tree number in a fChain
  
   /* Fixed size dimensions of array or collections stored in the TTree if any.
    * If this is too small, Segfaults may follow. */
   static const Int_t kMaxfPrtcls = 5000;
  
   //PrtclEvent        *event;
   Int_t           fPrtcls_;
   Double_t        fPrtcls_fP4_fCoordinates_fX[kMaxfPrtcls];   //[fPrtcls_]
   Double_t        fPrtcls_fP4_fCoordinates_fY[kMaxfPrtcls];   //[fPrtcls_]
   Double_t        fPrtcls_fP4_fCoordinates_fZ[kMaxfPrtcls];   //[fPrtcls_]
   Double_t        fPrtcls_fP4_fCoordinates_fT[kMaxfPrtcls];   //[fPrtcls_]
  
   Int_t           fPrtcls_fPDGCode[kMaxfPrtcls];   //[fPrtcls_]
   Int_t           fPrtcls_fChargeTimes3[kMaxfPrtcls];   //[fPrtcls_]
  
   Int_t          fPrtcls_fAnalysisStatus[kMaxfPrtcls];   //[fPrtcls_]

   // List of branches
   TBranch        *b_event_fPrtcls_;   //!
   TBranch        *b_fPrtcls_fP4_fCoordinates_fX;   //!
   TBranch        *b_fPrtcls_fP4_fCoordinates_fY;   //!
   TBranch        *b_fPrtcls_fP4_fCoordinates_fZ;   //!
   TBranch        *b_fPrtcls_fP4_fCoordinates_fT;   //!
  
   TBranch        *b_fPrtcls_fPDGCode;   //!
   TBranch        *b_fPrtcls_fChargeTimes3;   //!
   
   TBranch        *b_fPrtcls_fAnalysisStatus;   //!

//////////
// Output:
//////////
   TFile          *fOutFile;
   TTree          *fOutTree;
   TBranch        *fJetBranch;
   JetEvent       *fjEvent;

//////////
// Others:
//////////
   Timer mTimer;

   /* Energy counters: */ 
   TLorentzVector mPiPlus, mPiMinus, mPi0Gamma, mGamma, 
      mKaPlus, mKaMinus, mKSZero, mKLZero, 
      mProton, mAproton, mNeutron, mAneutron,
      mLambda0, mSigma, mElec, mMuon,
      mOthers, mEtSum;
  
   vector<size_t> mFlavorIndices; 
  
   /* Weights etc.: */
   double mPartSum, mChargSum, mChargWSum, mChargW2Sum, mW2;
   int mQuarkJetCharge, mFlavour;
  
   /* Temporary fractions: */
   double mChf;
   double mNhf;
   double mPhf;
   double mElf;
   double mMuf;
  
   /* Temporary masses: */
   double mChm;
   double mNhm;
   double mPhm;
   double mElm;
   double mMum;
};

#endif // JETANALYSIS_H