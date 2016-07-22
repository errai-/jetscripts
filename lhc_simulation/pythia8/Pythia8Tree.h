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
// Last modification: 27.10.2015                                 //
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
#include <map>
#include <algorithm>
#include <stdexcept>
#include <iterator>

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
using std::map;
using std::cout;
using std::endl;
using std::cerr;
using std::runtime_error;

/* A selector class to speed up lepton+ selection in ttbar */
class TTBarSelector : public UserHooks
{
public:
    TTBarSelector() {}
    virtual ~TTBarSelector() {}

    /* In UserHooks this returns false */
    virtual bool canVetoProcessLevel() { return true; }

    /* Returns true, if the event is vetoed (i.e. amount of hard proc leptons != 2) */
    virtual bool doVetoProcessLevel(Event& process);
};

class Pythia8Tree
{
public:
    /* The event loop should be called once after a proper initialization */
    void EventLoop();

    /* No copying or constructing with another instance */
    Pythia8Tree(const Pythia8Tree& other) = delete;
    Pythia8Tree& operator=(const Pythia8Tree&) = delete;

    /* Run settings are provided through the initializer */
    Pythia8Tree(string settings, string fileName, int mode);
    Pythia8Tree() :
                    mEvent(mPythia.event),
                    mInitialized(false)
    {
        cerr << "Pythia8Tree is intended to be used only with the non-default initializer" << endl;
    }
    /* ROOT has an awful behaviour with pointers. It should do the cleaning-up
       and such stuff well, but there is a ton of memory leaks that result
       already from ROOT 'being there'. At least the software runs - but the
       memory leaks and ROOT style pointer handling are regrettable. */
    ~Pythia8Tree() {}

protected:

    /* A handle for adding particle information */
    void                ParticleAdd(unsigned prt, int status);
    /* Particles needed by the hadronic flavor definition */
    void                GhostHadronAdd(unsigned prt);

    /* Loop over particles within an event: returns true if event is to be saved */
    bool                ParticleLoop();
    /* The logic within particleloop. */
    bool                ProcessParticle(unsigned prt);
    /* A custom phase within ProcessParticle */
    virtual int         CustomProcess(unsigned prt) { return 0; }

    /* See: HadronAndPartonSelector.cc in CMSSW. Indicates whether a ghost hadron
     * is in an excited state or not. Checks whether a hadron has a daughter of
     * the same flavour. Parameter quarkId is a PDG quark flavour. */
    bool                IsExcitedHadronState(unsigned prt, int quarkID);
    /* A function that checks whether a photon is originated from a pi0 and that
     * the energy of the photon-pair corresponds to the pion. returns 0 if
     * the origin is not a pion with good energy and 1 if it is */
    bool                GammaChecker(unsigned prt);
    /* Has the particle already been appended */
    bool                Absent(unsigned prt);

    /* Print the decay tree */
    void                PrintParticle(unsigned prt);
    void                PrintEvent();

    /* Corrected parton momentum */
    TLorentzVector      LastParton(unsigned prt);

protected:

    /* Indicator that the event loop can be run */
    bool                mInitialized;
    /* Whether or not store strange ghosts */
    bool                mUseStrange;
    /* A general-purpose counters for physics debugging */
    unsigned            mCounter;
    unsigned            mIterCount;

    unsigned            mHardProcCount;

    Pythia              mPythia;
    Event&              mEvent;

    TFile              *mFile;
    TTree              *mTree;
    TBranch            *mBranch;
    PrtclEvent         *mPrtclEvent;

    TTBarSelector       mTTBarSelector;

    int                 mNumErrs;
    int                 mNumEvents;
    int                 mMode;
    int                 mTimerStep;
    Timer               mTimer;

    vector<unsigned>    mSpecialIndices;
};


class P8GenericTree : public Pythia8Tree
{
public:
    P8GenericTree(string settings, string fileName, int mode) :
        Pythia8Tree(settings, fileName, mode) {}
    ~P8GenericTree() {}

protected:
    virtual int CustomProcess(unsigned prt);
};


class P8DijetTree : public Pythia8Tree
{
public:
    P8DijetTree(string settings, string fileName, int mode) :
        Pythia8Tree(settings, fileName, mode) {}
    ~P8DijetTree() {}

protected:
    /* Dijet specific particle logic */
    virtual int CustomProcess(unsigned prt);
};


class P8GammajetTree : public Pythia8Tree
{
public:
    P8GammajetTree(string settings, string fileName, int mode) :
        Pythia8Tree(settings, fileName, mode) {}
    ~P8GammajetTree() {}

protected:
    /* A handle for adding a hard process photon descended from prt */
    bool GammaAdd(unsigned prt);
    /* Gammajet specific particle logic */
    virtual int CustomProcess(unsigned prt);
};


class P8ZmumujetTree : public Pythia8Tree
{
public:
    P8ZmumujetTree(string settings, string fileName, int mode) :
        Pythia8Tree(settings, fileName, mode) {}
    ~P8ZmumujetTree() {}

protected:
    /* A handle for adding the two muons originating from a hard process Z prt */
    bool MuonAdd(unsigned prt);
    /* Zmumujet specific particle logic */
    virtual int CustomProcess(unsigned prt);
};


class P8ttbarjetTree : public Pythia8Tree
{
public:
    P8ttbarjetTree(string settings, string fileName, int mode) :
        Pythia8Tree(settings, fileName, mode) {}
    ~P8ttbarjetTree() {}

protected:
    /* A handle for adding the produced leptons in ttbar events */
    bool LeptonAdd(unsigned prt);
    /* ttbarjet specific particle logic */
    virtual int CustomProcess(unsigned prt);
};


#endif // PYTHIA8TREE
