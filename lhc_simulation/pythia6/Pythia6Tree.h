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
// Last modification: 19.11.2015                                 //
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
#include <map>
#include <algorithm>
#include <stdexcept>

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
using std::cerr;
using std::runtime_error;
using std::stoi;
using std::map;

class Pythia6Tree
{
public:
    /* The event loop should be called once after a proper initialization */
    void EventLoop();

    /* No copying or constructing with another instance */
    Pythia6Tree( const Pythia6Tree& other ) = delete;
    Pythia6Tree& operator=( const Pythia6Tree& ) = delete;

    /* Run settings are provided via the initializer */
    Pythia6Tree(Int_t nEvent, string fileName, Int_t nameId, const int mode);
    Pythia6Tree() :
        mInitialized(false)
    {
        cerr << "Pythia8Tree is intended to be used only with the non-default initializer" << endl;
    }
    /* ROOT has an awful behaviour with pointers. It should do the cleaning-up
       and such stuff well, but there is a ton of memory leaks that result
       already from ROOT 'being there'. At least the software runs - but the
       memory leaks and ROOT style pointer handling are regrettable. */
    ~Pythia6Tree() {}

protected:

    /* A handle for adding particle information */
    void                            ParticleAdd(unsigned prt,int status);
    /* Particles needed by the hadronic flavor definition */
    void                            GhostHadronAdd(unsigned prt, bool useStrange = false);

    /* Loop over particles within an event: return true if event is to be saved */
    bool                            ParticleLoop();
    /* The logic within particleloop */
    bool                            ProcessParticle(unsigned prt);

    /* See: HadronAndPartonSelector.cc in CMSSW. Indicates whether a ghost hadron
     * is in an excited state or not. Checks whether a hadron has a daughter of
     * the same flavour. Parameter quarkId is a PDG quark flavour. */
    bool                            IsExcitedHadronState(unsigned prt, int quarkId);
    /* A function that checks whether a photon is originated from a pi0 and that
     * the energy of the photon-pair corresponds to the pion. returns 0 if
     * the origin is not a pion with good energy and 1 if it is */
    bool                            GammaChecker(unsigned);
    /* A function for indicating the flavour history */
    void                            PropagateHistory(unsigned prt, int hard_prt);
    /* Has the particle already been appended */
    bool                            Absent(unsigned prt);
    /* A shortcut for finding out what stuff is partons (no beam remnants) */
    bool                            IsParton(unsigned prt);
    /* Check whether this is an FS parton */
    bool                            IsFSParton(unsigned prt);

    /* Settings that depend on the selected mode */
    void                            ModeSettings();
    /* General settings that are always used */
    void                            GeneralSettings();

    /* A handle for adding a hard process photon descended from the signal event photon */
    bool                            GammaAdd();
    /* A handle for adding the two muons originating from a hard process Z prt */
    bool                            MuonAdd();
    /* A handle for adding the produced leptons in ttbar events */
    bool                            LeptonAdd(unsigned prt);

    /* Print the decay tree */
    void                            PrintParticle(unsigned prt);
    void                            PrintEvent();

    /* Corrected parton momentum */
    TLorentzVector                  LastParton(unsigned prt);

protected:

    /* Indicator that the event loop can be run */
    bool                            mInitialized;
    /* A general-purpose counter for physics debugging */
    unsigned                        mCounter;

    TPythia6                       *mPythia;

    TFile                          *mFile;
    TTree                          *mTree;
    TBranch                        *mBranch;
    PrtclEvent                     *mPrtclEvent;

    int                             mNumEvents;
    int                             mMode;
    int                             mTimerStep;
    Timer                           mTimer;

    vector<unsigned>                mSpecialIndices;
    map<unsigned,TLorentzVector>    mPartonHistory;

    int                             mEnergy;
    int                             mTune;
    bool                            mMPI;
    bool                            mISR;
    bool                            mFSR;
    bool                            mUseStrange;

    static const unsigned           mNumSeeds = 40;
    const unsigned                  mSeeds[mNumSeeds] = {840744607,
                                                         431166825,
                                                         11489507,
                                                         859341684,
                                                         719632152,
                                                         384411333,
                                                         90405435,
                                                         297596781,
                                                         620424940,
                                                         829585206,
                                                         350220548,
                                                         862060943,
                                                         865146589,
                                                         11119376,
                                                         706126850,
                                                         761335296,
                                                         286390445,
                                                         408256820,
                                                         447625541,
                                                         368022699,
                                                         281922559,
                                                         852542479,
                                                         509348179,
                                                         175162098,
                                                         688006297,
                                                         512118632,
                                                         676751467,
                                                         212155085,
                                                         158795947,
                                                         68988051,
                                                         258456879,
                                                         625579469,
                                                         146828216,
                                                         582720998,
                                                         226158642,
                                                         439232438,
                                                         366169042,
                                                         745702146,
                                                         412672564,
                                                         177882235};
};

#endif // PYTHIA6TREE
