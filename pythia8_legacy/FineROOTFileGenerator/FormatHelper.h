#ifndef FORMATHELPER
#define FORMATHELPER

#include <TTree.h>
#include <TROOT.h>
#include <TObject.h>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

class EvtHdr {
  public:
    Int_t mIsPVgood;
    Int_t mHCALNoise;
    Long64_t mRun;
    Long64_t mEvent;
    Long64_t mLumi;
    Int_t mBunch;
    Int_t mNVtx;
    Int_t mNVtxGood;
    Int_t mOOTPUEarly;
    Int_t mOOTPULate;
    Int_t mINTPU;
    Int_t mNBX;
    Int_t mPVndof;
    Int_t mTrPu;
    Double_t mPVx;
    Double_t mPVy;
    Double_t mPVz;
    Double_t mBSx;
    Double_t mBSy;
    Double_t mBSz;
    Double_t mPthat;
    Double_t mWeight;
    Double_t mCaloRho;
    Double_t mPFRho;
};

class CaloMet {
  public:
    Double_t et_;
    Double_t sumEt_;
    Double_t phi_;
};

class PFMet {
  public:
    Double_t et_;
    Double_t sumEt_;
    Double_t phi_;
};

class Coordinates {
  public:
    Int_t jets;
    Double_t *fX; // [jets]
    Double_t *fY; // [jets]
    Double_t *fZ; // [jets]
    Double_t *fT; // [jets]

    Coordinates() {
      fX = 0; fY = 0; fZ = 0; fT = 0;
    }

    ~Coordinates() {
      delete [] fX;
      delete [] fY;
      delete [] fZ;
      delete [] fT;
    }

    void Refresh(Int_t jets_) {
      jets = jets_;
      if (fX) { delete [] fX; delete [] fY; delete [] fZ; delete [] fT;}
      if (!jets_) { fX = 0; fY = 0; fZ = 0; fT = 0; return; }
      fX = new Double_t[jets];
      fY = new Double_t[jets];
      fZ = new Double_t[jets];
      fT = new Double_t[jets];
    }
};

class GenJets {
  public:
    Coordinates fCoordinates;
};

class P4 {
  public:
    Coordinates fCoordinates;
};

class Caardinates {
  public:
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fT;

    Caardinates() {
      fX = 0; fY = 0; fZ = 0; fT = 0;
    }

    ~Caardinates() {}
};

class P3 {
  public:
    Caardinates fCoordinates;
    P3(){ }
    ~P3(){ }
};

class CaloJets {
  public:
    Int_t jets;
    P4 P4_;
    P4 genP4_;
    Double_t *genR_; // [jets]
    Double_t *cor_; // [jets]
    Double_t *unc_; // [jets]
    Double_t *uncSrc_; // [jets]
    Double_t *area_; // [jets]
    Double_t *looseID_; // [jets]
    Double_t *tightID_; // [jets]
    Double_t *emf_; // [jets]
    Double_t *fHPD_; // [jets]
    Double_t *fRBX_; // [jets]
    Double_t *n90hits_; // [jets]
    Double_t *nTrkCalo_; // [jets]
    Double_t *nTrkVtx_; // [jets]

    CaloJets() {
      jets = 0; tightID_ = 0;
    }

    ~CaloJets() {
      delete [] tightID_;
    }
};

class PFJets {
  public:
    P4 P4_;
    P4 genP4_;
    Double_t genR_;
    Double_t cor_;
    Double_t unc_;
    Double_t uncSrc_;
    Double_t area_;
    Double_t looseID_;
    Double_t tightID_;
    Double_t chf_;
    Double_t nhf_; 
    Double_t phf_; 
    Double_t elf_;  
    Double_t muf_;  
    Double_t hf_hf_; 
    Double_t hf_phf_; 
    Double_t hf_hm_;
    Double_t hf_phm_;
    Double_t chm_;
    Double_t nhm_;
    Double_t phm_;
    Double_t elm_;
    Double_t mum_;
    Double_t ncand_;
    Double_t beta_;
    Double_t betaStar_;

    PFJets() {
      tightID_ = 0; chf_ = 0; nhf_ = 0; phf_ = 0;
      elf_ = 0; muf_ = 0; hf_hf_ = 0; hf_phf_ = 0;
    }

    ~PFJets() {}
};

class FatJets {
  public:
    Int_t jets;
    P4 P4_;
    P4 genP4_;
    Double_t *genR_; // [jets]
    Double_t *cor_; // [jets]
    Double_t *unc_; // [jets]
    Double_t *uncSrc_; // [jets]
    Double_t *area_; // [jets]
    Double_t *looseID_; // [jets]
    Double_t *tightID_; // [jets]

    FatJets() {
      jets = 0; tightID_ = 0;    
    }

    ~FatJets() {
      delete [] tightID_;
    }
};

class FormatHelper {
  public:
    FormatHelper(){ }
    ~FormatHelper(){ }

    EvtHdr EvtHdr_;
    CaloMet CaloMet_;
    PFMet PFMet_;
    Int_t TriggerDecision_;
    Double_t L1Prescale_;
    Double_t HLTPrescale_;
    GenJets GenJets_;
    CaloJets CaloJets_;
    std::vector<PFJets> PFJets_;
    FatJets FatJets_;
    double PFJets__;
};

#endif
