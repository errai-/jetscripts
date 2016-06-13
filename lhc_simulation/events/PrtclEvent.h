//////////////////////////////////////////////////////////////////////
//                                                                  //
// A class structure for storing particle data produced with event  //
// generators (e.g. Pythia6/Pythia8/Herwig++) into ROOT trees.      //
//                                                                  //
// IMPORTANT NOTES:                                                 //
//                                                                  //
//  - fAnalysisStatus is a purely user-defined index for indicating //
//    why a particle has been saved.                                //
//     1: final state particles                                     //
//     2: final state particles which have interesting properties   //
//     3: outgoing hard process partons                             //
//     4: partonic final state partons                              //
//     5: ghost s hadrons                                           //
//     6: ghost c hadrons                                           //
//     7: ghost b hadrons                                           //
//     8: outgoing hard process partons with corrected momenta      //
//                                                                  //
//  - This class structure includes some "bad c++" because of       //
//    the way how ROOT works. Thus for instance some missing        //
//    semicolons should not be corrected.                           //
//                                                                  //
// Author: Hannu Siikonen (errai- @GitHub)                          //
// Last modification: 21.8.2015                                     //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#ifndef PRTCLEVENT_H
#define PRTCLEVENT_H

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

    void SetPxPyPzE(double px,double py,double pz,double e);
    void SetParams(int id,int status,int history_flavor);
    
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
    int fHistoryFlavor;
   
public:
    ClassDef(PrtclData,1)
};


// An event class for the particle data that is stored
class PrtclEvent : public TObject 
{
public:
    PrtclEvent();
    virtual ~PrtclEvent() { Reset(); };

    void AddPrtcl(double px,double py,double pz,double e,int id,int stat, int hf);
    PrtclData *InitPrtcl();

    void Clear(Option_t *option ="");
    void Reset(Option_t *option ="");

    TClonesArray *GetParts() const { return fPrtcls; }
    double fWeight;

private:
    size_t fN_Prtcl; //! Present amount of particles in the tree 
    
    TClonesArray *fPrtcls;
    static TClonesArray *fgPrtcls;

public:
    ClassDef(PrtclEvent,1)
};


#endif // PRTCLEVENT_H
