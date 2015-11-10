/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 10.11.2015                   //
/////////////////////////////////////////////////

#ifndef HARDPROC_H
#define HARDPROC_H

#include "JetBase.h"

using std::cout;
using std::endl;
using std::vector;
using std::cerr;

class HardProc : public JetBase 
{
public : 
    HardProc(TTree * tree, const char* file1, const char* file2, Int_t mode, Int_t definition) :
        JetBase::JetBase(tree, file1, file2, mode, definition)
    {
        fAddNonJet = false;
    }

    virtual void EventProcessing(Long64_t);
    
protected:
    vector<fastjet::PseudoJet> fCutJetParts;
};

#endif // HARDPROC_H