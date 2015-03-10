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
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using std::vector;
using std::cout;
using std::endl;

class JetData : public TObject {
public:
  // Use a pure ROOT LorentzVector so that for instance Pt can be found out even 
  // without the sources of this event class. This is a slightly better format than
  // TLorentzVector and is in use for instance in the KKousouris scripts (indirectly,
  // through CMSSW).
  ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > fP;
 
  Double_t fChf;
  Double_t fNhf;
  Double_t fPhf;
  Double_t fElf;
  Double_t fMuf;
  
  Double_t fChm;
  Double_t fNhm;
  Double_t fPhm;
  Double_t fElm;
  Double_t fMum;
  
  Char_t fFlavour;
  
  JetData() { Class()->IgnoreTObjectStreamer(); }
  virtual ~JetData() { }
  
  Double_t P() const { return fP.P(); }
  Double_t Pt() const { return fP.Pt(); }
  Double_t Eta() const { return fP.Eta(); }
  Double_t Phi() const { return fP.Phi(); }
  Double_t Mass() const { return fP.M(); }  
  
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

  
void Build(double,double,double,double,double,double,double,double,double,
  double, double, double, double, double,char);
  void Clear(Option_t *option ="");
  void Reset(Option_t *option ="");

  void SetNjet(Int_t n) { fNjet = n; }
  Int_t GetNjet() { return fNjet; }
  
  JetData *AddJet();

  TClonesArray *GetJets() const {return fJets;}
  
  ClassDef(JetEvent, 1)
};


#endif // JETEVENT_H
