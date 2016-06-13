/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 10.11.2015                   //
/////////////////////////////////////////////////

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "JetBase.h"

using std::cout;
using std::endl;
using std::vector;
using std::cerr;

class Analysis : public JetBase 
{
public : 
    Analysis(TTree * tree, const char* file1, const char* file2, Int_t mode, Int_t definition) :
        JetBase::JetBase(tree, file1, file2, mode, definition)
    {
        fParticleStudy = true;
    }
    
    virtual Bool_t   IsHadron(Int_t);
    virtual Bool_t   IsCharged(Int_t);
    virtual void     TypeSort();
    
    virtual void     Cuts();
    virtual Double_t PTD();
    virtual Double_t Sigma2();

    virtual void     PostProcessing(unsigned = 0);
    
protected:
    vector<fastjet::PseudoJet> fCutJetParts;
};

#endif // ANALYSIS_H