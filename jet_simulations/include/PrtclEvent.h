#ifndef PARTICLEEVENT_H
#define PARTICLEEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing particle data from simulations. //
// Hannu Siikonen 7.3.2015                                           //
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



class PrtclData : public TObject {
public:
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fE;
  
  Int_t fPDGCode;
  Int_t fChargeTimes3;
  
  Bool_t IsPi0Photon;
  Bool_t IsJetFlavor;
  Bool_t IsExcitedState;
  
  PrtclData() { Class()->IgnoreTObjectStreamer(); }
  virtual ~PrtclData() { }
  
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
  
  ClassDef(PrtclData,1)
};


class ParticleEvent : public TObject {
private:
  
  Int_t fNpart; //! Not saved to a tree
  
  TClonesArray *fParts;
  
  static TClonesArray *fgParts;

public:
  
  ParticleEvent(size_t = 10000);
  virtual ~ParticleEvent();

  void Build(double,double,double,double,int,double,int=0,int=0,int=0);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNpart(Int_t n) { fNpart = n; }
  Int_t GetNpart() { return fNpart; }
  
  PrtclData *AddParticle();

  TClonesArray *GetParts() const {return fParts;}
  
  ClassDef(ParticleEvent, 1)
};


#endif // SIMEVENT_H
