#ifndef PRTCLEVENT_H
#define PRTCLEVENT_H

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
#include "TClass.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using std::vector;
using std::cout;
using std::endl;

class PrtclData : public TObject {
public:
  // Use a pure ROOT LorentzVector so that for instance Pt can be found out even 
  // without the sources of this event class. This is a slightly better format than
  // TLorentzVector and is in use for instance in the KKousouris scripts (indirectly,
  // through CMSSW).
  ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > fP4;
  
  Int_t fPDGCode;
  Int_t fChargeTimes3;
  
  Bool_t IsPi0Photon;
  Bool_t IsJetFlavor;
  Bool_t IsExcitedState;
  
  PrtclData() { Class()->IgnoreTObjectStreamer(); }
  virtual ~PrtclData() { }
  
  Double_t P() const { return fP4.P(); }
  Double_t Pt() const { return fP4.Pt(); }
  Double_t Eta() const { return fP4.Eta(); }
  Double_t Phi() const { return fP4.Phi(); }
  Double_t Mass() const { return fP4.M(); }
  
  ClassDef(PrtclData,1)
};

class PrtclEvent : public TObject {
private:
  Int_t fNpart; //! Not saved to a tree; present amount of particles in the tree
  Size_t fSizeLim; //! The maximal amount of particles within an event
  
  TClonesArray *fParts;
  
  static TClonesArray *fgParts;
public:
  PrtclEvent(size_t = 10000);
  virtual ~PrtclEvent() { Reset(); };

  void Build(double,double,double,double,int,double,int=0,int=0,int=0);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNpart(Int_t n) { fNpart = n; }
  Int_t GetNpart() { return fNpart; }
  
  PrtclData *AddParticle();

  TClonesArray *GetParts() const {return fParts;}
  
  ClassDef(PrtclEvent, 1)
};


#endif // PRTCLEVENT_H
