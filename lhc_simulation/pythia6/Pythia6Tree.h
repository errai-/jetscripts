///////////////////////////////////////////////////////////////////
//                                                               //
// Class structure for storing Pythia6 particle data into trees. //            
//                                                               //
// Modes of operation:                                           //
//                                                               //
//    0: Generic events                                          //
//       -Normal QCD production for general studies              //
//                                                               //
//    1: Dijet events                                            //
//       -Normal QCD production with dijet-specific settings     //
//                                                               //
//    2: Photon+Jet events                                       //
//       - Photon+Jet production                                 //
//                                                               //
//    3: Zmumu+Jet events                                        //
//       - Z+Jet production with a Z -> mumu setting             //
//                                                               //
//    4: Ttbarlepton+jet events                                  //
//       - Ttbar production with WW -> qqbarllbar                //
//                                                               //
// Author: Hannu Siikonen (errai- @GitHub)                       //
// Last modification: 21.8.2015                                  //
//                                                               //
///////////////////////////////////////////////////////////////////

#ifndef PYTHIA6TREE
#define PYTHIA6TREE 

/* Stdlib header file for input and output. */
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

/* Header file to access Pythia 8 program elements. */
#include "Pythia8/Pythia.h"

/* ROOT */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"

/* scripts */
#include "../generic/help_functions.h"
#include "../events/PrtclEvent.h"

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::runtime_error;


class Pythia6Tree
{
public:

    Pythia6Tree(Int_t nEvent, string fileName, Int_t nameId, const int mode) :
        mMode(mode), mNumEvents(nEvent)
    {
        /* Create an instance of the Pythia event generator: */
        mPythia = new TPythia6;
        /* Set a seed value according to the run index and make sure it is used: */
        mPythia->SetMRPY(1,10000*nameId);
        mPythia->SetMRPY(2,0);
        
        /* Event type: */
        ModeSettings();
        /* Other settings: */
        GeneralSettings();

        /* Try to create a file to write */
        mFile = TFile::Open(fileName.c_str(), "RECREATE");
        if(!mFile->IsOpen()) throw runtime_error("Creating an output file failed");
        mFile->SetCompressionLevel(1);

        /* Output tree: */
        mTree = new TTree("Pythia6Tree", "Pythia6 particle data.");
        if(!mTree) throw runtime_error("Creating a tree failed");
        mTree->SetAutoSave(100000000); /* 0.1 GBytes */
        mTree->SetCacheSize(10000000); /* 100 MBytes */
        TTree::SetBranchStyle(1); /* New branch style */

        /* Connect an event to the tree */
        mPrtclEvent = new PrtclEvent();
        if (!mPrtclEvent) throw runtime_error("Creating an event handle failed");
        mBranch = mTree->Branch("event", &mPrtclEvent, 32000,4);
        if (!mBranch) throw runtime_error("Associating the event handle with the tree failed");
        mBranch->SetAutoDelete(kFALSE);
        mTree->BranchRef();
        
        /* Setup a custom event timer */
        mTimerStep = 1000;
        mTimer.setParams(mNumEvents,mTimerStep);       
        mTimer.startTiming();
    }
    
    ~Pythia6Tree() 
    {
        delete mPrtclEvent; mPrtclEvent = 0;
        delete mPythia; mPythia = 0;
        mFile->Close();
    };
    
    void ModeSettings();
    void GeneralSettings();
    
    void EventLoop();
    
    void ParticleAdd(std::size_t,int);
    
    bool GammaChecker(unsigned);
    
    bool ParticleLoop();
    
    void GammaAdd();
    void MuonAdd();
    void LeptonAdd(std::size_t);
    
protected:
    
    TPythia6 *mPythia;
    
    TFile *mFile;
    TTree *mTree;
    TBranch *mBranch;
    PrtclEvent *mPrtclEvent;
    vector< pair<std::size_t,int> > mCandidates;
    std::size_t mNextCand;
    
    int mNumEvents;
    int mMode;
    int mTimerStep;
    Timer mTimer;
    
    vector<std::size_t> mSpecialIndices;
};

#endif // PYTHIA6TREE