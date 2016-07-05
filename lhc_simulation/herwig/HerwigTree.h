////////////////////////////////////////////////////////////////////
//                                                                //
// Class structure for storing Herwig++ particle data into trees. //
//                                                                //
// Modes of operation:                                            //
//                                                                //
//    0: Generic events                                           //
//       -Normal QCD production for general studies               //
//                                                                //
//    1: Dijet events                                             //
//       -Normal QCD production with dijet-specific settings      //
//                                                                //
//    2: Photon+Jet events                                        //
//       - Photon+Jet production                                  //
//                                                                //
//    3: Zmumu+Jet events                                         //
//       - Z+Jet production with a Z -> mumu setting              //
//                                                                //
//    4: Ttbarlepton+jet events                                   //
//       - Ttbar production with WW -> qqbarllbar                 //
//                                                                //
// Author: Hannu Siikonen (errai- @GitHub)                        //
// Last modification: 30.6.2016                                   //
//                                                                //
////////////////////////////////////////////////////////////////////


#ifndef HERWIGTREE_H
#define HERWIGTREE_H

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Config/Unitsystem.h"

#include "../generic/help_functions.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "../events/PrtclEvent.h"

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <stdexcept>

using std::cout;
using std::endl;
using std::runtime_error;


namespace Herwig {

using namespace ThePEG;

class HerwigTree: public AnalysisHandler {

public:

    HerwigTree() {}

    /** Analyze a given Event.
     * @param event pointer to the Event to be analyzed.
     * @param ieve the event number.
     * @param loop the number of times this event has been presented.
     * @param state nonzero if the event has been manipulated. */
    virtual void analyze(tEventPtr event, long ieve, int loop, int state);

    /** Standard Init function, called exactly once. */
    static void Init() { static ClassDocumentation<HerwigTree> documentation("The HerwigTree class is intended for saving gen-level data into tree structures."); }

protected:

    virtual IBPtr clone() const { return new_ptr(*this); }
    virtual IBPtr fullclone() const { return new_ptr(*this); }

    virtual void doinitrun();
    virtual void dofinish();

    /** @name Functions for the particle loop. */
    //@{
    bool gammaAdd(tPPtr);

    bool muonAdd(tPPtr);

    bool leptonAdd(tPPtr);
    //@}

    void print_parents(const tPPtr&);

    /** @name Help methods for the analysis */
    //@{
   /** ThePEG does not provide useful status codes and the status has to be studied manually.
     * This method is a mock-up of the CMSSW-way to calculate the status code. However, these
     * codes are not at the moment in use and the method is here only for reference */
    int getStatusCode(const tPPtr&) const;

    int gammaChecker(const tPPtr&);

    int isExcitedHadronState(const tPPtr&, int);

    void particleAdd(const tPPtr&, int);
    //@}

    /** @name Variables for the analysis */
    //@{

    /* The hardest event */
    tcEventBasePtr mHard;

    TFile *mFile;
    TTree *mTree;
    TBranch *mBranch;
    PrtclEvent *mPrtclEvent;

    int mNumEvents;
    int mMode;
    int mTimerStep;
    Timer mTimer;

    vector<unsigned> mSpecialIndices;
    //@}

private:

   /** The static object used to initialize the description of this class. */
    static NoPIOClassDescription<HerwigTree> initHerwigTree;

   /** The assignment operator is private and must never be called. */
    HerwigTree & operator=(const HerwigTree &);

};

/* Initialization, closing and administrative stuff: */

void HerwigTree::doinitrun()
{
    AnalysisHandler::doinitrun();
    try {
        string fileName = "particles_herwig";

        /* In a general multithread-case, generate a thread-unique root file name */
        fileName += "_";
        fileName += generator()->runName();
        fileName += ".root";

        size_t pos = fileName.find("jet_");
        string modeName = fileName.substr(17,pos-17);
        if (modeName=="generic") {
            mMode = 0;
        } else if (modeName=="di") {
            mMode = 1;
        } else if (modeName=="gamma") {
            mMode = 2;
        } else if (modeName=="Z") {
            mMode = 3;
        } else if (modeName=="ttbar") {
            mMode = 4;
        } else {
            throw runtime_error("Bad mode");
        }

        /* Setup a root file */
        mFile = new TFile (fileName.c_str(),"RECREATE");
        if (!mFile) throw runtime_error("Creating an output file failed");
        mFile->SetCompressionLevel(1);

        /* Setup a root tree */
        mTree = new TTree ("HerwigTree","Herwig++ particle data.");
        if (!mTree) throw runtime_error("Creating a tree failed");
        mTree->SetAutoSave(100000000);  /* 0.1 GBytes */
        mTree->SetCacheSize(10000000);  /* 10 MBytes */
        TTree::SetBranchStyle(1); /* new style */

        /* Connect an event handle with the tree */
        mPrtclEvent = new PrtclEvent;
        if (!mPrtclEvent) throw runtime_error("Creating an event handle failed");
        mBranch = mTree->Branch("event", &mPrtclEvent, 32000,4);
        if (!mBranch) throw runtime_error("Associating the event handle with the tree failed");
        mBranch->SetAutoDelete(kFALSE);
        mTree->BranchRef();

        /* Timing functions */
        mTimerStep = 1000;
        mTimer.setParams(mNumEvents,mTimerStep);
        mTimer.startTiming();

    } catch (std::exception& e) {
        cout << "An error occurred: " << e.what() << endl;
    }
}

void HerwigTree::dofinish()
{
    AnalysisHandler::dofinish();

    mTree->GetCurrentFile();
    mTree->AutoSave("Overwrite");
    mFile->Close();

    cout << "A tree has been written into a .root file" << endl;
}

NoPIOClassDescription<HerwigTree> HerwigTree::initHerwigTree;
// Definition of the static class description member.

}



#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HerwigTree. */
template <>
struct BaseClassTrait<Herwig::HerwigTree,1> {
    /** Typedef of the first base class of HerwigTree. */
    typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HerwigTree class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HerwigTree>
    : public ClassTraitsBase<Herwig::HerwigTree> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::HerwigTree"; }
    /** Return the name(s) of the shared library (or libraries) be loaded to get
    *  access to the HerwigTree class and any other class on which it depends
    *  (except the base class). */
    static string library() { return "libHerwigTree.so"; }
};

/** @endcond */

}

#endif /* HERWIGTREE_H */
