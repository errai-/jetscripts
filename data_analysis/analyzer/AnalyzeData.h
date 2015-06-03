
#ifndef AnalyzeData_h
#define AnalyzeData_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>

#include <vector>
#include <iostream>
#include <string>
#include <cassert>
#include <set>

#include "AnalyzeHelper.h"
#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrector.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc"
#include "tdrstyle_mod.C"

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxGenJets_ = 1;
const Int_t kMaxCaloJets_ = 200;
const Int_t kMaxPFJets_ = 200;
const Int_t kMaxFatJets_ = 2;

class AnalyzeData {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain 
    Long64_t        loopLimit;
    Int_t           isMC;
    Int_t           isDT;     // As a default behaviour, isDT is simply !isMC
    Int_t           testing;  // A flag for testing memory and time usage

    // Declaration of leaf types
 //QCDEvent        *events;
    Bool_t          EvtHdr__mIsPVgood;
    Bool_t          EvtHdr__mHCALNoise;
    Int_t           EvtHdr__mRun;
    Int_t           EvtHdr__mEvent;
    Int_t           EvtHdr__mLumi;
    Int_t           EvtHdr__mBunch;
    Int_t           EvtHdr__mNVtx;
    Int_t           EvtHdr__mNVtxGood;
    Int_t           EvtHdr__mOOTPUEarly;
    Int_t           EvtHdr__mOOTPULate;
    Int_t           EvtHdr__mINTPU;
    Int_t           EvtHdr__mNBX;
    Float_t         EvtHdr__mPVndof;
    Float_t         EvtHdr__mTrPu;
    Float_t         EvtHdr__mPVx;
    Float_t         EvtHdr__mPVy;
    Float_t         EvtHdr__mPVz;
    Float_t         EvtHdr__mBSx;
    Float_t         EvtHdr__mBSy;
    Float_t         EvtHdr__mBSz;
    Float_t         EvtHdr__mPthat;
    Float_t         EvtHdr__mWeight;
    Float_t         EvtHdr__mCaloRho;
    Float_t         EvtHdr__mPFRho;
    Float_t         CaloMet__et_;
    Float_t         CaloMet__sumEt_;
    Float_t         CaloMet__phi_;
    Float_t         PFMet__et_;
    Float_t         PFMet__sumEt_;
    Float_t         PFMet__phi_;
    vector<int>     TriggerDecision_;
    vector<int>     L1Prescale_;
    vector<int>     HLTPrescale_;
    
    Int_t           GenJets__;
    Double_t        GenJets__fCoordinates_fX[kMaxGenJets_];   //[GenJets__]
    Double_t        GenJets__fCoordinates_fY[kMaxGenJets_];   //[GenJets__]
    Double_t        GenJets__fCoordinates_fZ[kMaxGenJets_];   //[GenJets__]
    Double_t        GenJets__fCoordinates_fT[kMaxGenJets_];   //[GenJets__]
    Int_t           CaloJets__;
    Double_t        CaloJets__P4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__P4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__P4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__P4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__genP4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__genP4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__genP4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
    Double_t        CaloJets__genP4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__genR_[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__cor_[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__unc_[kMaxCaloJets_];   //[CaloJets__]
    vector<float>   CaloJets__uncSrc_[kMaxCaloJets_];
    Float_t         CaloJets__area_[kMaxCaloJets_];   //[CaloJets__]
    Bool_t          CaloJets__looseID_[kMaxCaloJets_];   //[CaloJets__]
    Bool_t          CaloJets__tightID_[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__emf_[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__fHPD_[kMaxCaloJets_];   //[CaloJets__]
    Float_t         CaloJets__fRBX_[kMaxCaloJets_];   //[CaloJets__]
    Int_t           CaloJets__n90hits_[kMaxCaloJets_];   //[CaloJets__]
    Int_t           CaloJets__nTrkCalo_[kMaxCaloJets_];   //[CaloJets__]
    Int_t           CaloJets__nTrkVtx_[kMaxCaloJets_];   //[CaloJets__]
    Int_t           PFJets__;
    Double_t        PFJets__P4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__genR_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__cor_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__unc_[kMaxPFJets_];   //[PFJets__]
    vector<float>   PFJets__uncSrc_[kMaxPFJets_];
    Float_t         PFJets__area_[kMaxPFJets_];   //[PFJets__]
    Bool_t          PFJets__looseID_[kMaxPFJets_];   //[PFJets__]
    Bool_t          PFJets__tightID_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__chf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__nhf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__phf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__elf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__muf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__hf_hf_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__hf_phf_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__hf_hm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__hf_phm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__chm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__nhm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__phm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__elm_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__mum_[kMaxPFJets_];   //[PFJets__]
    Int_t           PFJets__ncand_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__beta_[kMaxPFJets_];   //[PFJets__]
    Float_t         PFJets__betaStar_[kMaxPFJets_];   //[PFJets__]
    Int_t           FatJets__;
    Double_t        FatJets__P4__fCoordinates_fX[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__P4__fCoordinates_fY[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__P4__fCoordinates_fZ[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__P4__fCoordinates_fT[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__genP4__fCoordinates_fX[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__genP4__fCoordinates_fY[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__genP4__fCoordinates_fZ[kMaxFatJets_];   //[FatJets__]
    Double_t        FatJets__genP4__fCoordinates_fT[kMaxFatJets_];   //[FatJets__]
    Float_t         FatJets__genR_[kMaxFatJets_];   //[FatJets__]
    Float_t         FatJets__cor_[kMaxFatJets_];   //[FatJets__]
    Float_t         FatJets__unc_[kMaxFatJets_];   //[FatJets__]
    vector<float>   FatJets__uncSrc_[kMaxFatJets_];
    Float_t         FatJets__area_[kMaxFatJets_];   //[FatJets__]
    Bool_t          FatJets__looseID_[kMaxFatJets_];   //[FatJets__]
    Bool_t          FatJets__tightID_[kMaxFatJets_];   //[FatJets__]

    // List of branches
    TBranch        *b_events_EvtHdr__mIsPVgood;   //!
    TBranch        *b_events_EvtHdr__mHCALNoise;   //!
    TBranch        *b_events_EvtHdr__mRun;   //!
    TBranch        *b_events_EvtHdr__mEvent;   //!
    TBranch        *b_events_EvtHdr__mLumi;   //!
    TBranch        *b_events_EvtHdr__mBunch;   //!
    TBranch        *b_events_EvtHdr__mNVtx;   //!
    TBranch        *b_events_EvtHdr__mNVtxGood;   //!
    TBranch        *b_events_EvtHdr__mOOTPUEarly;   //!
    TBranch        *b_events_EvtHdr__mOOTPULate;   //!
    TBranch        *b_events_EvtHdr__mINTPU;   //!
    TBranch        *b_events_EvtHdr__mNBX;   //!
    TBranch        *b_events_EvtHdr__mPVndof;   //!
    TBranch        *b_events_EvtHdr__mTrPu;   //!
    TBranch        *b_events_EvtHdr__mPVx;   //!
    TBranch        *b_events_EvtHdr__mPVy;   //!
    TBranch        *b_events_EvtHdr__mPVz;   //!
    TBranch        *b_events_EvtHdr__mBSx;   //!
    TBranch        *b_events_EvtHdr__mBSy;   //!
    TBranch        *b_events_EvtHdr__mBSz;   //!
    TBranch        *b_events_EvtHdr__mPthat;   //!
    TBranch        *b_events_EvtHdr__mWeight;   //!
    TBranch        *b_events_EvtHdr__mCaloRho;   //!
    TBranch        *b_events_EvtHdr__mPFRho;   //!
    TBranch        *b_events_CaloMet__et_;   //!
    TBranch        *b_events_CaloMet__sumEt_;   //!
    TBranch        *b_events_CaloMet__phi_;   //!
    TBranch        *b_events_PFMet__et_;   //!
    TBranch        *b_events_PFMet__sumEt_;   //!
    TBranch        *b_events_PFMet__phi_;   //!
    TBranch        *b_events_TriggerDecision_;   //!
    TBranch        *b_events_L1Prescale_;   //!
    TBranch        *b_events_HLTPrescale_;   //!
    TBranch        *b_events_GenJets__;   //!
    TBranch        *b_GenJets__fCoordinates_fX;   //!
    TBranch        *b_GenJets__fCoordinates_fY;   //!
    TBranch        *b_GenJets__fCoordinates_fZ;   //!
    TBranch        *b_GenJets__fCoordinates_fT;   //!
    TBranch        *b_events_CaloJets__;   //!
    TBranch        *b_CaloJets__P4__fCoordinates_fX;   //!
    TBranch        *b_CaloJets__P4__fCoordinates_fY;   //!
    TBranch        *b_CaloJets__P4__fCoordinates_fZ;   //!
    TBranch        *b_CaloJets__P4__fCoordinates_fT;   //!
    TBranch        *b_CaloJets__genP4__fCoordinates_fX;   //!
    TBranch        *b_CaloJets__genP4__fCoordinates_fY;   //!
    TBranch        *b_CaloJets__genP4__fCoordinates_fZ;   //!
    TBranch        *b_CaloJets__genP4__fCoordinates_fT;   //!
    TBranch        *b_CaloJets__genR_;   //!
    TBranch        *b_CaloJets__cor_;   //!
    TBranch        *b_CaloJets__unc_;   //!
    TBranch        *b_CaloJets__uncSrc_;   //!
    TBranch        *b_CaloJets__area_;   //!
    TBranch        *b_CaloJets__looseID_;   //!
    TBranch        *b_CaloJets__tightID_;   //!
    TBranch        *b_CaloJets__emf_;   //!
    TBranch        *b_CaloJets__fHPD_;   //!
    TBranch        *b_CaloJets__fRBX_;   //!
    TBranch        *b_CaloJets__n90hits_;   //!
    TBranch        *b_CaloJets__nTrkCalo_;   //!
    TBranch        *b_CaloJets__nTrkVtx_;   //!
    TBranch        *b_events_PFJets__;   //!
    TBranch        *b_PFJets__P4__fCoordinates_fX;   //!
    TBranch        *b_PFJets__P4__fCoordinates_fY;   //!
    TBranch        *b_PFJets__P4__fCoordinates_fZ;   //!
    TBranch        *b_PFJets__P4__fCoordinates_fT;   //!
    TBranch        *b_PFJets__genP4__fCoordinates_fX;   //!
    TBranch        *b_PFJets__genP4__fCoordinates_fY;   //!
    TBranch        *b_PFJets__genP4__fCoordinates_fZ;   //!
    TBranch        *b_PFJets__genP4__fCoordinates_fT;   //!
    TBranch        *b_PFJets__genR_;   //!
    TBranch        *b_PFJets__cor_;   //!
    TBranch        *b_PFJets__unc_;   //!
    TBranch        *b_PFJets__uncSrc_;   //!
    TBranch        *b_PFJets__area_;   //!
    TBranch        *b_PFJets__looseID_;   //!
    TBranch        *b_PFJets__tightID_;   //!
    TBranch        *b_PFJets__chf_;   //!
    TBranch        *b_PFJets__nhf_;   //!
    TBranch        *b_PFJets__phf_;   //!
    TBranch        *b_PFJets__elf_;   //!
    TBranch        *b_PFJets__muf_;   //!
    TBranch        *b_PFJets__hf_hf_;   //!
    TBranch        *b_PFJets__hf_phf_;   //!
    TBranch        *b_PFJets__hf_hm_;   //!
    TBranch        *b_PFJets__hf_phm_;   //!
    TBranch        *b_PFJets__chm_;   //!
    TBranch        *b_PFJets__nhm_;   //!
    TBranch        *b_PFJets__phm_;   //!
    TBranch        *b_PFJets__elm_;   //!
    TBranch        *b_PFJets__mum_;   //!
    TBranch        *b_PFJets__ncand_;   //!
    TBranch        *b_PFJets__beta_;   //!
    TBranch        *b_PFJets__betaStar_;   //!
    TBranch        *b_events_FatJets__;   //!
    TBranch        *b_FatJets__P4__fCoordinates_fX;   //!
    TBranch        *b_FatJets__P4__fCoordinates_fY;   //!
    TBranch        *b_FatJets__P4__fCoordinates_fZ;   //!
    TBranch        *b_FatJets__P4__fCoordinates_fT;   //!
    TBranch        *b_FatJets__genP4__fCoordinates_fX;   //!
    TBranch        *b_FatJets__genP4__fCoordinates_fY;   //!
    TBranch        *b_FatJets__genP4__fCoordinates_fZ;   //!
    TBranch        *b_FatJets__genP4__fCoordinates_fT;   //!
    TBranch        *b_FatJets__genR_;   //!
    TBranch        *b_FatJets__cor_;   //!
    TBranch        *b_FatJets__unc_;   //!
    TBranch        *b_FatJets__uncSrc_;   //!
    TBranch        *b_FatJets__area_;   //!
    TBranch        *b_FatJets__looseID_;   //!
    TBranch        *b_FatJets__tightID_;   //!

    AnalyzeData(TTree * = 0, Long64_t = 0, Int_t = 0, Int_t = 0);
    virtual ~AnalyzeData();
    virtual Int_t    GetEntry(Long64_t);
    virtual Long64_t LoadTree(Long64_t);
    virtual void     Init(TTree *);
    virtual void     Loop(string);
    virtual void     Show(Long64_t = -1);
};

#endif

#ifdef AnalyzeData_cxx

AnalyzeData::AnalyzeData(TTree *tree, Long64_t events, Int_t mc, Int_t tests) : fChain(0) 
{
    loopLimit=events;
    isMC=mc;
    isDT=!isMC;
    testing=tests;
    Init(tree);
}

AnalyzeData::~AnalyzeData()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t AnalyzeData::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t AnalyzeData::LoadTree(Long64_t entry)
{
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

void AnalyzeData::Init(TTree *tree)
{
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("EvtHdr_.mIsPVgood", &EvtHdr__mIsPVgood, &b_events_EvtHdr__mIsPVgood);
    fChain->SetBranchAddress("EvtHdr_.mHCALNoise", &EvtHdr__mHCALNoise, &b_events_EvtHdr__mHCALNoise);
    fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
    fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
    fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
    fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch, &b_events_EvtHdr__mBunch);
    fChain->SetBranchAddress("EvtHdr_.mNVtx", &EvtHdr__mNVtx, &b_events_EvtHdr__mNVtx);
    fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood, &b_events_EvtHdr__mNVtxGood);
    fChain->SetBranchAddress("EvtHdr_.mOOTPUEarly", &EvtHdr__mOOTPUEarly, &b_events_EvtHdr__mOOTPUEarly);
    fChain->SetBranchAddress("EvtHdr_.mOOTPULate", &EvtHdr__mOOTPULate, &b_events_EvtHdr__mOOTPULate);
    fChain->SetBranchAddress("EvtHdr_.mINTPU", &EvtHdr__mINTPU, &b_events_EvtHdr__mINTPU);
    fChain->SetBranchAddress("EvtHdr_.mNBX", &EvtHdr__mNBX, &b_events_EvtHdr__mNBX);
    fChain->SetBranchAddress("EvtHdr_.mPVndof", &EvtHdr__mPVndof, &b_events_EvtHdr__mPVndof);
    fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu, &b_events_EvtHdr__mTrPu);
    fChain->SetBranchAddress("EvtHdr_.mPVx", &EvtHdr__mPVx, &b_events_EvtHdr__mPVx);
    fChain->SetBranchAddress("EvtHdr_.mPVy", &EvtHdr__mPVy, &b_events_EvtHdr__mPVy);
    fChain->SetBranchAddress("EvtHdr_.mPVz", &EvtHdr__mPVz, &b_events_EvtHdr__mPVz);
    fChain->SetBranchAddress("EvtHdr_.mBSx", &EvtHdr__mBSx, &b_events_EvtHdr__mBSx);
    fChain->SetBranchAddress("EvtHdr_.mBSy", &EvtHdr__mBSy, &b_events_EvtHdr__mBSy);
    fChain->SetBranchAddress("EvtHdr_.mBSz", &EvtHdr__mBSz, &b_events_EvtHdr__mBSz);
    fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat, &b_events_EvtHdr__mPthat);
    fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);
    fChain->SetBranchAddress("EvtHdr_.mCaloRho", &EvtHdr__mCaloRho, &b_events_EvtHdr__mCaloRho);
    fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho, &b_events_EvtHdr__mPFRho);
    fChain->SetBranchAddress("CaloMet_.et_", &CaloMet__et_, &b_events_CaloMet__et_);
    fChain->SetBranchAddress("CaloMet_.sumEt_", &CaloMet__sumEt_, &b_events_CaloMet__sumEt_);
    fChain->SetBranchAddress("CaloMet_.phi_", &CaloMet__phi_, &b_events_CaloMet__phi_);
    fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
    fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);
    fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_, &b_events_PFMet__phi_);
    fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
    fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
    fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
    fChain->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
    fChain->SetBranchAddress("GenJets_.fCoordinates.fX", &GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
    fChain->SetBranchAddress("GenJets_.fCoordinates.fY", &GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
    fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", &GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
    fChain->SetBranchAddress("GenJets_.fCoordinates.fT", &GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
    fChain->SetBranchAddress("CaloJets_", &CaloJets__, &b_events_CaloJets__);
    fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fX", CaloJets__P4__fCoordinates_fX, &b_CaloJets__P4__fCoordinates_fX);
    fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fY", CaloJets__P4__fCoordinates_fY, &b_CaloJets__P4__fCoordinates_fY);
    fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fZ", CaloJets__P4__fCoordinates_fZ, &b_CaloJets__P4__fCoordinates_fZ);
    fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fT", CaloJets__P4__fCoordinates_fT, &b_CaloJets__P4__fCoordinates_fT);
    fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fX", CaloJets__genP4__fCoordinates_fX, &b_CaloJets__genP4__fCoordinates_fX);
    fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fY", CaloJets__genP4__fCoordinates_fY, &b_CaloJets__genP4__fCoordinates_fY);
    fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fZ", CaloJets__genP4__fCoordinates_fZ, &b_CaloJets__genP4__fCoordinates_fZ);
    fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fT", CaloJets__genP4__fCoordinates_fT, &b_CaloJets__genP4__fCoordinates_fT);
    fChain->SetBranchAddress("CaloJets_.genR_", CaloJets__genR_, &b_CaloJets__genR_);
    fChain->SetBranchAddress("CaloJets_.cor_", CaloJets__cor_, &b_CaloJets__cor_);
    fChain->SetBranchAddress("CaloJets_.unc_", CaloJets__unc_, &b_CaloJets__unc_);
    fChain->SetBranchAddress("CaloJets_.uncSrc_", CaloJets__uncSrc_, &b_CaloJets__uncSrc_);
    fChain->SetBranchAddress("CaloJets_.area_", CaloJets__area_, &b_CaloJets__area_);
    fChain->SetBranchAddress("CaloJets_.looseID_", CaloJets__looseID_, &b_CaloJets__looseID_);
    fChain->SetBranchAddress("CaloJets_.tightID_", CaloJets__tightID_, &b_CaloJets__tightID_);
    fChain->SetBranchAddress("CaloJets_.emf_", CaloJets__emf_, &b_CaloJets__emf_);
    fChain->SetBranchAddress("CaloJets_.fHPD_", CaloJets__fHPD_, &b_CaloJets__fHPD_);
    fChain->SetBranchAddress("CaloJets_.fRBX_", CaloJets__fRBX_, &b_CaloJets__fRBX_);
    fChain->SetBranchAddress("CaloJets_.n90hits_", CaloJets__n90hits_, &b_CaloJets__n90hits_);
    fChain->SetBranchAddress("CaloJets_.nTrkCalo_", CaloJets__nTrkCalo_, &b_CaloJets__nTrkCalo_);
    fChain->SetBranchAddress("CaloJets_.nTrkVtx_", CaloJets__nTrkVtx_, &b_CaloJets__nTrkVtx_);
    fChain->SetBranchAddress("PFJets_", &PFJets__, &b_events_PFJets__);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fX", PFJets__P4__fCoordinates_fX, &b_PFJets__P4__fCoordinates_fX);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fY", PFJets__P4__fCoordinates_fY, &b_PFJets__P4__fCoordinates_fY);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fZ", PFJets__P4__fCoordinates_fZ, &b_PFJets__P4__fCoordinates_fZ);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fT", PFJets__P4__fCoordinates_fT, &b_PFJets__P4__fCoordinates_fT);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fX", PFJets__genP4__fCoordinates_fX, &b_PFJets__genP4__fCoordinates_fX);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fY", PFJets__genP4__fCoordinates_fY, &b_PFJets__genP4__fCoordinates_fY);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fZ", PFJets__genP4__fCoordinates_fZ, &b_PFJets__genP4__fCoordinates_fZ);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fT", PFJets__genP4__fCoordinates_fT, &b_PFJets__genP4__fCoordinates_fT);
    fChain->SetBranchAddress("PFJets_.genR_", PFJets__genR_, &b_PFJets__genR_);
    fChain->SetBranchAddress("PFJets_.cor_", PFJets__cor_, &b_PFJets__cor_);
    fChain->SetBranchAddress("PFJets_.unc_", PFJets__unc_, &b_PFJets__unc_);
    fChain->SetBranchAddress("PFJets_.uncSrc_", PFJets__uncSrc_, &b_PFJets__uncSrc_);
    fChain->SetBranchAddress("PFJets_.area_", PFJets__area_, &b_PFJets__area_);
    fChain->SetBranchAddress("PFJets_.looseID_", PFJets__looseID_, &b_PFJets__looseID_);
    fChain->SetBranchAddress("PFJets_.tightID_", PFJets__tightID_, &b_PFJets__tightID_);
    fChain->SetBranchAddress("PFJets_.chf_", PFJets__chf_, &b_PFJets__chf_);
    fChain->SetBranchAddress("PFJets_.nhf_", PFJets__nhf_, &b_PFJets__nhf_);
    fChain->SetBranchAddress("PFJets_.phf_", PFJets__phf_, &b_PFJets__phf_);
    fChain->SetBranchAddress("PFJets_.elf_", PFJets__elf_, &b_PFJets__elf_);
    fChain->SetBranchAddress("PFJets_.muf_", PFJets__muf_, &b_PFJets__muf_);
    fChain->SetBranchAddress("PFJets_.hf_hf_", PFJets__hf_hf_, &b_PFJets__hf_hf_);
    fChain->SetBranchAddress("PFJets_.hf_phf_", PFJets__hf_phf_, &b_PFJets__hf_phf_);
    fChain->SetBranchAddress("PFJets_.hf_hm_", PFJets__hf_hm_, &b_PFJets__hf_hm_);
    fChain->SetBranchAddress("PFJets_.hf_phm_", PFJets__hf_phm_, &b_PFJets__hf_phm_);
    fChain->SetBranchAddress("PFJets_.chm_", PFJets__chm_, &b_PFJets__chm_);
    fChain->SetBranchAddress("PFJets_.nhm_", PFJets__nhm_, &b_PFJets__nhm_);
    fChain->SetBranchAddress("PFJets_.phm_", PFJets__phm_, &b_PFJets__phm_);
    fChain->SetBranchAddress("PFJets_.elm_", PFJets__elm_, &b_PFJets__elm_);
    fChain->SetBranchAddress("PFJets_.mum_", PFJets__mum_, &b_PFJets__mum_);
    fChain->SetBranchAddress("PFJets_.ncand_", PFJets__ncand_, &b_PFJets__ncand_);
    fChain->SetBranchAddress("PFJets_.beta_", PFJets__beta_, &b_PFJets__beta_);
    fChain->SetBranchAddress("PFJets_.betaStar_", PFJets__betaStar_, &b_PFJets__betaStar_);
    fChain->SetBranchAddress("FatJets_", &FatJets__, &b_events_FatJets__);
    fChain->SetBranchAddress("FatJets_.P4_.fCoordinates.fX", FatJets__P4__fCoordinates_fX, &b_FatJets__P4__fCoordinates_fX);
    fChain->SetBranchAddress("FatJets_.P4_.fCoordinates.fY", FatJets__P4__fCoordinates_fY, &b_FatJets__P4__fCoordinates_fY);
    fChain->SetBranchAddress("FatJets_.P4_.fCoordinates.fZ", FatJets__P4__fCoordinates_fZ, &b_FatJets__P4__fCoordinates_fZ);
    fChain->SetBranchAddress("FatJets_.P4_.fCoordinates.fT", FatJets__P4__fCoordinates_fT, &b_FatJets__P4__fCoordinates_fT);
    fChain->SetBranchAddress("FatJets_.genP4_.fCoordinates.fX", FatJets__genP4__fCoordinates_fX, &b_FatJets__genP4__fCoordinates_fX);
    fChain->SetBranchAddress("FatJets_.genP4_.fCoordinates.fY", FatJets__genP4__fCoordinates_fY, &b_FatJets__genP4__fCoordinates_fY);
    fChain->SetBranchAddress("FatJets_.genP4_.fCoordinates.fZ", FatJets__genP4__fCoordinates_fZ, &b_FatJets__genP4__fCoordinates_fZ);
    fChain->SetBranchAddress("FatJets_.genP4_.fCoordinates.fT", FatJets__genP4__fCoordinates_fT, &b_FatJets__genP4__fCoordinates_fT);
    fChain->SetBranchAddress("FatJets_.genR_", FatJets__genR_, &b_FatJets__genR_);
    fChain->SetBranchAddress("FatJets_.cor_", FatJets__cor_, &b_FatJets__cor_);
    fChain->SetBranchAddress("FatJets_.unc_", FatJets__unc_, &b_FatJets__unc_);
    fChain->SetBranchAddress("FatJets_.uncSrc_", FatJets__uncSrc_, &b_FatJets__uncSrc_);
    fChain->SetBranchAddress("FatJets_.area_", FatJets__area_, &b_FatJets__area_);
    fChain->SetBranchAddress("FatJets_.looseID_", FatJets__looseID_, &b_FatJets__looseID_);
    fChain->SetBranchAddress("FatJets_.tightID_", FatJets__tightID_, &b_FatJets__tightID_);
}

void AnalyzeData::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}

#endif // #ifdef AnalyzeData_cxx
