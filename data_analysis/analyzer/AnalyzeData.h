
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

const Int_t kMaxPFJets_ = 200;

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
    Int_t           EvtHdr__mRun;
    Int_t           EvtHdr__mEvent;
    Int_t           EvtHdr__mLumi;
    Int_t           EvtHdr__mNVtxGood;
    Float_t         EvtHdr__mTrPu;
    Float_t         EvtHdr__mPthat;
    Float_t         EvtHdr__mPFRho;
    Float_t         PFMet__et_;
    Float_t         PFMet__sumEt_;
    vector<int>     TriggerDecision_;
    vector<int>     L1Prescale_;
    vector<int>     HLTPrescale_;
    Int_t           PFJets__;
    Double_t        PFJets__P4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__P4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
    Double_t        PFJets__genP4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
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
    Float_t         PFJets__betaStar_[kMaxPFJets_];   //[PFJets__]

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

    fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun);
    fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent);
    fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi);
    fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood);
    fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu);
    fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat);
    fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho);
    fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_);
    fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_);
    fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_);
    fChain->SetBranchAddress("L1Prescale_", &L1Prescale_);
    fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_);
    fChain->SetBranchAddress("PFJets_", &PFJets__);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fX", PFJets__P4__fCoordinates_fX);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fY", PFJets__P4__fCoordinates_fY);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fZ", PFJets__P4__fCoordinates_fZ);
    fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fT", PFJets__P4__fCoordinates_fT);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fX", PFJets__genP4__fCoordinates_fX);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fY", PFJets__genP4__fCoordinates_fY);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fZ", PFJets__genP4__fCoordinates_fZ);
    fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fT", PFJets__genP4__fCoordinates_fT);
    fChain->SetBranchAddress("PFJets_.area_", PFJets__area_);
    fChain->SetBranchAddress("PFJets_.looseID_", PFJets__looseID_);
    fChain->SetBranchAddress("PFJets_.tightID_", PFJets__tightID_);
    fChain->SetBranchAddress("PFJets_.chf_", PFJets__chf_);
    fChain->SetBranchAddress("PFJets_.nhf_", PFJets__nhf_);
    fChain->SetBranchAddress("PFJets_.phf_", PFJets__phf_);
    fChain->SetBranchAddress("PFJets_.elf_", PFJets__elf_);
    fChain->SetBranchAddress("PFJets_.muf_", PFJets__muf_);
    fChain->SetBranchAddress("PFJets_.hf_hf_", PFJets__hf_hf_);
    fChain->SetBranchAddress("PFJets_.hf_phf_", PFJets__hf_phf_);
    fChain->SetBranchAddress("PFJets_.betaStar_", PFJets__betaStar_);
    
    /* First disable all branches and then enable the branches in use */
    fChain->SetBranchStatus("*",0);
    fChain->SetBranchStatus("EvtHdr_.mRun", 1);
    fChain->SetBranchStatus("EvtHdr_.mEvent", 1);
    fChain->SetBranchStatus("EvtHdr_.mLumi", 1);
    fChain->SetBranchStatus("EvtHdr_.mNVtxGood", 1);
    fChain->SetBranchStatus("EvtHdr_.mTrPu", 1);
    fChain->SetBranchStatus("EvtHdr_.mPthat", 1);
    fChain->SetBranchStatus("EvtHdr_.mPFRho", 1);
    fChain->SetBranchStatus("PFMet_.et_", 1);
    fChain->SetBranchStatus("PFMet_.sumEt_", 1);
    fChain->SetBranchStatus("TriggerDecision_", 1);
    fChain->SetBranchStatus("L1Prescale_", 1);
    fChain->SetBranchStatus("HLTPrescale_", 1);
    fChain->SetBranchStatus("PFJets_", 1);
    fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fX",1);
    fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fY",1);
    fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fZ",1);
    fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fT",1);
    fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fX", 1);
    fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fY", 1);
    fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fZ", 1);
    fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fT", 1);
    fChain->SetBranchStatus("PFJets_.area_", 1);
    fChain->SetBranchStatus("PFJets_.looseID_", 1);
    fChain->SetBranchStatus("PFJets_.tightID_", 1);
    fChain->SetBranchStatus("PFJets_.chf_", 1);
    fChain->SetBranchStatus("PFJets_.nhf_", 1);
    fChain->SetBranchStatus("PFJets_.phf_", 1);
    fChain->SetBranchStatus("PFJets_.elf_", 1);
    fChain->SetBranchStatus("PFJets_.muf_", 1);
    fChain->SetBranchStatus("PFJets_.hf_hf_", 1);
    fChain->SetBranchStatus("PFJets_.hf_phf_", 1);
    fChain->SetBranchStatus("PFJets_.betaStar_", 1);

    
}

void AnalyzeData::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}

#endif // #ifdef AnalyzeData_cxx
