#ifndef SIMEVENT_H
#define SIMEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing particle data from simulations. //
// Hannu Siikonen 13.2.2015                                          //
// (special thanks to Rene Brun's ROOT examples)                     //
///////////////////////////////////////////////////////////////////////

//#include "TROOT.h"
#include <vector>
#include <iostream>
#include <fstream>
// #include <cstdint>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TArrayD.h"
#include "TClonesArray.h"
#include "TProcessID.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TClass.h"

using std::vector;
using std::cout;
using std::endl;



class SimParticle : public TObject {
public:
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fE;
  
  Int_t fPDGCode;
  
  Bool_t IsPi0Photon;
  Bool_t IsJetFlavor;
  Bool_t IsExcitedState;
  
  SimParticle() { Class()->IgnoreTObjectStreamer(); }
  virtual ~SimParticle() { }
  
  Double_t P() const { return pow( pow(fPx,2) + pow(fPy,2) + pow(fPz,2), 0.5); }
  Double_t Pt() const { return pow( pow(fPx,2) + pow(fPy,2), 0.5 ); }
  Double_t Eta() const {
    if ( P() - fPz == 0 ) return 1000000000000;
    else return 0.5*TMath::Log( ( P() + fPz )/( P() - fPz ) ); 
  }
  Double_t Phi() const { return TMath::ATan2( fPy, fPx ); }
  Double_t Mass() const { return pow( pow(fE,2) - pow( P(), 2 ), 0.5 ); }
  
  TLorentzVector GetLorentz() const {
    TLorentzVector tmpVect(fPx,fPy,fPz,fE);
    return tmpVect;
  }
  
  ClassDef(SimParticle,1)
};


class SimEvent : public TObject {
private:
  
  Int_t fNpart; //! Not saved to a tree
  
  TClonesArray *fParts;
  
  static TClonesArray *fgParts;

public:
  
  SimEvent(size_t = 10000);
  virtual ~SimEvent();

  void Build(double,double,double,double,int,int=0,int=0,int=0);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNpart(Int_t n) { fNpart = n; }
  Int_t GetNpart() { return fNpart; }
  
  SimParticle *AddParticle();

  TClonesArray *GetParts() const {return fParts;}
  
  ClassDef(SimEvent, 1)
};


#endif // SIMEVENT_H
