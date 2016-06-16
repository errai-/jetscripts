///////////////////////////////////////////////////////////////////////
// A generic event class for storing jet data from simulations.      //
// Hannu Siikonen 19.5.2015                                          //
///////////////////////////////////////////////////////////////////////

#ifndef JETEVENT_H
#define JETEVENT_H

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

struct JetVariables
{
    float chf;
    float nhf;
    float phf;
    float elf;
    float muf;

    float chm;
    float nhm;
    float phm;
    float elm;
    float mum;

    //float nextDR;
    float partonPT;
    float matchPT;

    int constituents;
    float PTD;
    float Sigma2;

    float DR;
    float Alpha;
    float DPhi;

    void SetZero() {
        chf = 0;
        nhf = 0;
        phf = 0;
        elf = 0;
        muf = 0;
        chm = 0;
        nhm = 0;
        phm = 0;
        elm = 0;
        mum = 0;
        partonPT = 0;
        matchPT = 0;
        constituents = 0;
        PTD = 0;
        Sigma2 = 0;
        DR = 0;
        Alpha = 0;
        DPhi = 0;
        //nextDR = 0;
    }
};


// The information concerning one jet, stored to a ROOT tree.
class JetData : public TObject {
public:
    JetData() { }
    virtual ~JetData() { }

    void SetPxPyPzE(double,double,double,double);
    void SetParams(JetVariables&,int);

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

    int   fFlav;

    float fChf;
    float fNhf;
    float fPhf;
    float fElf;
    float fMuf;

    //float fChm;
    //float fNhm;
    //float fPhm;
    //float fElm;
    //float fMum;

    float fPartonPT;
    float fMatchPT;

    int   fConstituents;
    float fPTD;
    float fSigma2;

    float fDR;
    //float fDR_Next;
    float fAlpha;
    float fDPhi;

public:
   ClassDef(JetData,1)
};


// An event class for the jet data that is stored
class JetEvent : public TObject {
public:
    JetEvent(size_t = 10000);
    virtual ~JetEvent();

    void AddJet(double,double,double,double,JetVariables&,double,int);
    JetData *InitJet();

    void Clear(Option_t *option ="");
    void Reset(Option_t *option ="");

private:
    size_t fN_Jet; //!
    size_t fSizeLim; //! Maximal amount of particles within an event

    TClonesArray *fJets;
    static TClonesArray *fgJets;
    double fWeight;

public:
    ClassDef(JetEvent,1)
};


#endif // JETEVENT_H
