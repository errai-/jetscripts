#ifndef PRTCLEVENT_H
#define PRTCLEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing particle data from simulations. //
// Hannu Siikonen 19.5.2015                                          //
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


// The information concerning one particle, stored to a ROOT tree.
class PrtclData : public TObject 
{
public:
    
    PrtclData() { }
    virtual ~PrtclData() { }

    void SetPxPyPzE(double,double,double,double);
    void SetParams(int,int);
    
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
    
    int fPDGCode;
    int fAnalysisStatus;
   
public:
    ClassDef(PrtclData,1)
};


// An event class for the particle data that is stored
class PrtclEvent : public TObject 
{
public:
    PrtclEvent(size_t = 10000);
    virtual ~PrtclEvent() { Reset(); };

    void AddPrtcl(double,double,double,double,int,int);
    PrtclData *InitPrtcl();

    void Clear(Option_t *option ="");
    void Reset(Option_t *option ="");

    TClonesArray *GetParts() const { return fPrtcls; }
    double fWeight;

private:
    size_t fN_Prtcl; //! Present amount of particles in the tree 
    size_t fSizeLim; //! Maximal amount of particles within an event
    
    TClonesArray *fPrtcls;
    static TClonesArray *fgPrtcls;

public:
    ClassDef(PrtclEvent,1)
};


#endif // PRTCLEVENT_H
