///////////////////////////////////////////////////////////////////
//                                                               //
// Class structure for storing Pythia8 particle data into trees. //            
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

#ifndef PYTHIA8TREE
#define PYTHIA8TREE 

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

/* scripts */
#include "../generic/help_functions.h"
#include "../events/PrtclEvent.h"

using namespace Pythia8;
using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::runtime_error;

class Pythia8Tree
{
public:

    Pythia8Tree(string settings, string fileName, int mode) : mEvent(mPythia.event)
    {
        /* Initialization of the Pythia8 run */
        if (!mPythia.readFile(settings.c_str())) throw std::invalid_argument("Error while reading settings"); 
        if (!mPythia.init()) throw runtime_error("Pythia8 initialization failed");
        mPythia.settings.listChanged();

        mNumEvents = mPythia.mode("Main:numberOfEvents");
        mMode = mode;

        /* Try to create a file to write */
        mFile = new TFile(fileName.c_str(), "RECREATE");
        if(!mFile->IsOpen()) throw runtime_error("Creating an output file failed");
        mFile->SetCompressionLevel(1);

        /* Create a tree. Autosave every 100 Mb, cache of 10 Mb */
        mTree = new TTree("Pythia8Tree","Pythia8 particle data.");
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
    
    ~Pythia8Tree() 
    {
        delete mPrtclEvent; mPrtclEvent = 0;
        mFile->Close();
    };
    
    void EventLoop();
    
    void ParticleAdd(std::size_t,int);
    TLorentzVector Vogel(unsigned);
    
    bool ParticleLoop();
    
    bool GammaAdd(std::size_t);
    bool MuonAdd(std::size_t);
    bool LeptonAdd(std::size_t);
    
    bool GammaChecker(std::size_t);
    
    bool IsExcitedHadronState(std::size_t, int);
    void GhostHadronAdd(std::size_t, bool = false);
    
protected:
    
    Pythia mPythia;
    Event& mEvent;
    
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

#endif // PYTHIA8TREE