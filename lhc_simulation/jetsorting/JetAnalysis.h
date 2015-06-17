/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 12.3.2015                    //
/////////////////////////////////////////////////

#ifndef JETANALYSIS_H
#define JETANALYSIS_H

// General
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector2.h>
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
#include "../events/PrtclEvent.h"
#include "../events/JetEvent.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include "../generic/help_functions.h"
#include "../generic/tdrstyle_mod1.C"

using std::cout;
using std::endl;
using std::vector;

class JetAnalysis 
{
public : 
    JetAnalysis(TTree *, const char*, const char*, Int_t, Int_t);

    ~JetAnalysis();
    
    virtual void     Init(TTree*); // Chain
    virtual void     InitCI();
    virtual void     InitFP();
    
    virtual Int_t    GetEntry(Long64_t);
    virtual Long64_t LoadTree(Long64_t);
    virtual void     Show(Long64_t = -1);
    
    virtual void     EventLoop();
    virtual void     ParticlesToJetsorterInput();
    virtual void     JetLoop();
    virtual Bool_t   GoodEvent();
    virtual Bool_t   IsHadron(Int_t);
    virtual Bool_t   IsCharged(Int_t);
    virtual void     PhysicsFlavor(size_t);
    virtual void     FlavorLoop(size_t);
    virtual Int_t    ChargeSign(Int_t);
    virtual void     ParticleLoop(size_t);
    virtual void     TypeSort();
    virtual void     FillerHandle( vector<TProfile*> &, Double_t, Double_t );
    virtual void     HistFill(int);    
    
    virtual void     Cuts();
    virtual Double_t PTD();
    virtual Double_t Sigma2();
    
    virtual void WriteResults();

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
    vector<fastjet::PseudoJet> sortedJets, jetParts, cutJetParts, 
                               fjInputs, hiddenInputs;
////////
// Root:
////////

  static const int ptBins = 48.;
  const double ptRange[ptBins+1]=
  //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
  //3637, 3832, 4037};//  
/////////
// Input:
/////////  
    TTree          *fChain;   //! 
    Int_t           fCurrent; //! current Tree number in a fChain
    
    /* Fixed size dimensions of array or collections stored in the TTree if any.
     * If this is too small, Segfaults may follow. */
    static const Int_t kMaxfPrtcls = 5000;
    
    //PrtclEvent        *event;
    Double_t        fWeight;
    Int_t           fPrtcls_;
    Double_t        fX[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fY[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fZ[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fT[kMaxfPrtcls];   //[fPrtcls_]
    
    Int_t           fPDGCode[kMaxfPrtcls];   //[fPrtcls_]
    Int_t           fChargeTimes3[kMaxfPrtcls];   //[fPrtcls_]
    
    Int_t           fAnalysisStatus[kMaxfPrtcls];   //[fPrtcls_]

    // List of branches
    TBranch        *b_fWeight;   //!
    TBranch        *b_fPrtcls_;   //!
    TBranch        *b_fX;   //!
    TBranch        *b_fY;   //!
    TBranch        *b_fZ;   //!
    TBranch        *b_fT;   //!

    TBranch        *b_fPDGCode;   //!
    TBranch        *b_fChargeTimes3;   //!
    
    TBranch        *b_fAnalysisStatus;   //!

//////////
// Output:
//////////
    TFile          *fOutFile;
    TFile          *fOutFile2;
    string         fOutFileName;
    TTree          *fOutTree;
    TBranch        *fJetBranch;
    JetEvent       *fjEvent;

//////////
// Graphs:
//////////
    TProfile* gluonQuark;
    vector<TH1D*> chargeIndicator;
    
    // Book histograms.
    vector<TProfile*> fractionProfilesGluon;
    vector<TProfile*> fractionProfilesQuark;
    vector<TProfile*> fractionProfilesLQuark;
    vector<TProfile*> fractionProfilesHQuark;
    vector<TProfile*> fractionProfilesAll;
    
//////////
// Others:
//////////
    Timer mTimer;

    /* Energy counters: */ 
    TLorentzVector mPiPlus, mPiMinus, mPi0Gamma, mGamma, 
        mKaPlus, mKaMinus, mKSZero, mKLZero, 
        mProton, mAproton, mNeutron, mAneutron,
        mLambda0, mSigma, mElec, mMuon, mOthers, mEtSum;
    
    vector<size_t> mFlavorIndices; 
    
    Int_t           mMode;       /* Event type */
    Int_t           mDefinition; /* Flavour definition */
    
    vector<Int_t>   mMuonList; /* Z-jets */
    Int_t           mGammaId;    /* gamma-jets */
    vector<Int_t>   mPartonList; /* Physics definition */
 
    /* Weights etc.: */
    double mChargSum, mChargWSum, mChargW2Sum, mW2;
    int mQuarkJetCharge, mFlavour, mIsHadron;
    
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