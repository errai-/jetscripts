#ifndef SIMEVENT_H
#define SIMEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing jet data from simulations.      //
// Hannu Siikonen 23.2.2015                                          //
// (special thanks to Rene Brun's ROOT examples)                     //
///////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cassert>

#include "TROOT.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TProcessID.h"
#include "TMath.h"
#include "TClass.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


// The information concerning one jet, stored to a ROOT tree.
class JetData : public TObject {
public:
   JetData() { Class()->IgnoreTObjectStreamer(); }
   virtual ~JetData() { }
   
   void SetPxPyPzE(double,double,double,double);
   void SetParams(double,double,double,double,double,double,double,double,
      double,double,int);
   
   Double_t P() const { return fP4.P(); }
   Double_t Pt() const { return fP4.Pt(); }
   Double_t Eta() const { return fP4.Eta(); }
   Double_t Phi() const { return fP4.Phi(); }
   Double_t Mass() const { return fP4.M(); }  

private:
   /* Use a pure ROOT LorentzVector so that for instance Pt can be found out
    * even without the sources of this event class. This is a slightly better 
    * format than TLorentzVector and is in use for instance in the KKousouris 
    * scripts (indirectly, through CMSSW). */ 
   ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > fP4;
   
   double fChf;
   double fNhf;
   double fPhf;
   double fElf;
   double fMuf;
   
   double fChm;
   double fNhm;
   double fPhm;
   double fElm;
   double fMum;
   
   int fFlav;
 
public:
   ClassDef(JetData,1)
};


// An event class for the jet data that is stored
class JetEvent : public TObject {
public:
   JetEvent(size_t = 1000);
   virtual ~JetEvent();

   void AddJet(double,double,double,double,double,double,double,double,double,
      double, double, double, double, double,int);
   JetData *InitJet();
   
   void Clear(Option_t *option ="");
   void Reset(Option_t *option ="");

private:
   size_t fN_Jet; //! 
   size_t fSizeLim; //! Maximal amount of particles within an event

   TClonesArray *fJets;
   static TClonesArray *fgJets;
 
public:
   ClassDef(JetEvent, 1)
};


#endif // JETEVENT_H
