/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 10.11.2015                   //
/////////////////////////////////////////////////

#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include "JetBase.h"

using std::cout;
using std::endl;
using std::vector;
using std::cerr;

class Histograms : public JetBase 
{
public : 
    Histograms(TTree * tree, const char* file1, const char* file2, Int_t mode, Int_t definition) :
        JetBase::JetBase(tree, file1, file2, mode, definition)
    {
        fParticleStudy = true;
        fJetCuts = false;
        fOpposite = false;
        gluonQuark = new TProfile("gq","gq",ptBins,ptRange);  
        InitFP();
    }
    
    virtual void     InitFP();
    virtual void     Finalize();

    virtual void     PostProcessing(unsigned = 0);
        
    virtual Int_t    ChargeSign(Int_t);
    virtual void     FillerHandle( vector<TProfile*> &, Double_t, Double_t );
    virtual void     HistFill(int);
    
    virtual void     WriteResults();
    
protected:
    vector<PseudoJet> fCutJetParts;
    int fQuarkJetCharge;
    bool fOpposite;
    
//////////
// Graphs:
//////////
    
  static const int ptBins = 48.;
  const double ptRange[ptBins+1]=
  //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
  //3637, 3832, 4037};//  
    
    TProfile* gluonQuark;
    
    // Book histograms.
    vector<TProfile*> fractionProfilesGluon;
    vector<TProfile*> fractionProfilesQuark;
    vector<TProfile*> fractionProfilesLQuark;
    vector<TProfile*> fractionProfilesHQuark;
    vector<TProfile*> fractionProfilesAll;
};

#endif // HISTOGRAMS_H