/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 10.11.2015                   //
/////////////////////////////////////////////////

#ifndef TEST_H
#define TEST_H

#include "JetBase.h"

#include <utility>

using std::cout;
using std::endl;
using std::vector;
using std::cerr;
using std::pair;

class Test : public JetBase 
{
public : 
    Test(TTree * tree, const char* file1, const char* file2, Int_t mode, Int_t definition) :
        JetBase::JetBase(tree, file1, file2, mode, definition)
    {
        fAccepted = 0;
        fAddNonJet = true;
        fNoisy = false;
    }

    virtual bool JetLoop();
    virtual void Finalize();
    
protected:
    vector<PseudoJet> fCutJetParts;
    int fAccepted;
    bool fNoisy;
};

#endif // TEST_H