#ifndef SIMEVENT_H
#define SIMEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing jet data from simulations.      //
// Hannu Siikonen 23.2.2015                                          //
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

class JetData : public TObject {
public:
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fE;
 
  Double_t fChf;
  Double_t fNhf;
  Double_t fPhf;
  Double_t fElf;
  Double_t fMuf;
  
  Char_t fFlavour;
  
  JetData() { Class()->IgnoreTObjectStreamer(); }
  virtual ~JetData() { }
  
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
  
  ClassDef(JetData,1)
};


class JetEvent : public TObject {
private:
  
  Int_t fNjet; //! Not saved to a tree
  
  TClonesArray *fJets;
  
  static TClonesArray *fgJets;

public:
  
  JetEvent(size_t = 10000);
  virtual ~JetEvent();

  void Build(double,double,double,double,double,double,double,double,double,char);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNjet(Int_t n) { fNjet = n; }
  Int_t GetNjet() { return fNjet; }
  
  JetData *AddJet();

  TClonesArray *GetParts() const {return fJets;}
  
  ClassDef(JetEvent, 1)
};


#endif // JETEVENT_H
