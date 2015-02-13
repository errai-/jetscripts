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
  
  Double_t Pt() const;
  Double_t Eta() const;
  Double_t Phi() const;
  
  SimParticle() { }
  virtual ~SimParticle() { }
  
  ClassDef(SimParticle,1);
};

class SimEvent : public TObject {
private:
  
  Int_t fNpart;
  
  TClonesArray *fParts;
  
  static TClonesArray *fgParts;

public:
  
  SimEvent(size_t = 10000);
  virtual ~SimEvent();

  void Build(double,double,double,double,int);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNpart(Int_t n) { fNpart = n; }
  Int_t GetNpart() const { return fNpart; }
  
  SimParticle *AddParticle();

  TClonesArray *GetParts() const {return fParts;}
  
  ClassDef(SimEvent, 1);
};


#endif // SIMEVENT_H
