
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
    Int_t           mRun;
    Int_t           mEvent;
    Int_t           mLumi;
    Int_t           mNVtxGood;
    Float_t         mTrPu;
    Float_t         mPthat;
    Float_t         mPFRho;
    Float_t         PFMet__et_;
    Float_t         PFMet__sumEt_;
    vector<int>     TriggerDecision_;
    vector<int>     L1_;
    vector<int>     HLT_;
    Int_t           PFJets_;
    Double_t        fWeight;
    Double_t        fX[kMaxPFJets_];   //[PFJets_]
    Double_t        fY[kMaxPFJets_];   //[PFJets_]
    Double_t        fZ[kMaxPFJets_];   //[PFJets_]
    Double_t        fT[kMaxPFJets_];   //[PFJets_]
    Double_t        gen_fX[kMaxPFJets_];   //[PFJets_]
    Double_t        gen_fY[kMaxPFJets_];   //[PFJets_]
    Double_t        gen_fZ[kMaxPFJets_];   //[PFJets_]
    Double_t        gen_fT[kMaxPFJets_];   //[PFJets_]
    Float_t         area_[kMaxPFJets_];   //[PFJets_]
    Bool_t          looseID_[kMaxPFJets_];   //[PFJets_]
    Bool_t          tightID_[kMaxPFJets_];   //[PFJets_]
    Float_t         chf_[kMaxPFJets_];   //[PFJets_]
    Float_t         nhf_[kMaxPFJets_];   //[PFJets_]
    Float_t         phf_[kMaxPFJets_];   //[PFJets_]
    Float_t         elf_[kMaxPFJets_];   //[PFJets_]
    Float_t         muf_[kMaxPFJets_];   //[PFJets_]
    Double_t         chf[kMaxPFJets_];   //[PFJets_]
    Double_t         nhf[kMaxPFJets_];   //[PFJets_]
    Double_t         phf[kMaxPFJets_];   //[PFJets_]
    Double_t         elf[kMaxPFJets_];   //[PFJets_]
    Double_t         muf[kMaxPFJets_];   //[PFJets_]
    Float_t         hf_hf_[kMaxPFJets_];   //[PFJets_]
    Float_t         hf_phf_[kMaxPFJets_];   //[PFJets_]
    Float_t         betaStar_[kMaxPFJets_];   //[PFJets_]

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
    
    if (isMC!=3) {
        fChain->SetBranchAddress("EvtHdr_.mRun", &mRun);
        fChain->SetBranchAddress("EvtHdr_.mEvent", &mEvent);
        fChain->SetBranchAddress("EvtHdr_.mLumi", &mLumi);
        fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &mNVtxGood);
        fChain->SetBranchAddress("EvtHdr_.mTrPu", &mTrPu);
        fChain->SetBranchAddress("EvtHdr_.mPthat", &mPthat);
        fChain->SetBranchAddress("EvtHdr_.mPFRho", &mPFRho);
        fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_);
        fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_);
        fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_);
        fChain->SetBranchAddress("L1Prescale_", &L1_);
        fChain->SetBranchAddress("HLTPrescale_", &HLT_);
        fChain->SetBranchAddress("PFJets_", &PFJets_);
        fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fX", fX);
        fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fY", fY);
        fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fZ", fZ);
        fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fT", fT);
        fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fX", gen_fX);
        fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fY", gen_fY);
        fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fZ", gen_fZ);
        fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fT", gen_fT);
        fChain->SetBranchAddress("PFJets_.area_", area_);
        fChain->SetBranchAddress("PFJets_.looseID_", looseID_);
        fChain->SetBranchAddress("PFJets_.tightID_", tightID_);
        fChain->SetBranchAddress("PFJets_.chf_", chf_);
        fChain->SetBranchAddress("PFJets_.nhf_", nhf_);
        fChain->SetBranchAddress("PFJets_.phf_", phf_);
        fChain->SetBranchAddress("PFJets_.elf_", elf_);
        fChain->SetBranchAddress("PFJets_.muf_", muf_);
        fChain->SetBranchAddress("PFJets_.hf_hf_", hf_hf_);
        fChain->SetBranchAddress("PFJets_.hf_phf_", hf_phf_);
        fChain->SetBranchAddress("PFJets_.betaStar_", betaStar_);
    } else {
        fChain->SetBranchAddress("fJets", &PFJets_);
        fChain->SetBranchAddress("fWeight", &fWeight);
        fChain->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
        fChain->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
        fChain->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
        fChain->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
        fChain->SetBranchAddress("fJets.fChf", chf);
        fChain->SetBranchAddress("fJets.fNhf", nhf);
        fChain->SetBranchAddress("fJets.fPhf", phf);
        fChain->SetBranchAddress("fJets.fElf", elf);
        fChain->SetBranchAddress("fJets.fMuf", muf);
    }
    
    /* First disable all branches and then enable the branches in use */
    fChain->SetBranchStatus("*",0);
    if (isMC!=3) {
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
    } else {
        fChain->SetBranchStatus("fJets",1);
        fChain->SetBranchStatus("fWeight",1);
        fChain->SetBranchStatus("fJets.fP4.fCoordinates.fX",1);
        fChain->SetBranchStatus("fJets.fP4.fCoordinates.fY",1);
        fChain->SetBranchStatus("fJets.fP4.fCoordinates.fZ",1);
        fChain->SetBranchStatus("fJets.fP4.fCoordinates.fT",1);
        fChain->SetBranchStatus("fJets.fChf",1);
        fChain->SetBranchStatus("fJets.fNhf",1);
        fChain->SetBranchStatus("fJets.fPhf",1);
        fChain->SetBranchStatus("fJets.fElf",1);
        fChain->SetBranchStatus("fJets.fMuf",1);
    }
    
}

void AnalyzeData::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}

#endif // #ifdef AnalyzeData_cxx
