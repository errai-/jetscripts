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
// Last modification: 20.8.2015                                   //
//                                                                //
////////////////////////////////////////////////////////////////////

#ifndef STOREPARTICLES_H
#define STOREPARTICLES_H

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Config/Unitsystem.h"

#include "../generic/help_functions.h"

#include "TTree.h"
#include "TFile.h"
#include "../events/PrtclEvent.h"

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cassert>
using std::cout;
using std::endl;

namespace jetanalysis {

using namespace ThePEG;

class HerwigppTree: public AnalysisHandler {

public:

    HerwigppTree() {}
    
    ~HerwigppTree() {}
    
    /** Analyze a given Event.
        * @param event pointer to the Event to be analyzed.
        * @param ieve the event number.
        * @param loop the number of times this event has been presented.
        * @param state nonzero if the event has been manipulated. */
    virtual void analyze(tEventPtr event, long ieve, int loop, int state);

    /** Standard Init function, called exactly once. */
    static void Init() { static ClassDocumentation<HerwigppTree> documentation("No documentation"); }

protected:
   
    /** @name Generic dummy-methods. */
    //@{
    IBPtr clone() const { return new_ptr(*this); }
    IBPtr fullclone() const { return new_ptr(*this); }
    //@}

    /** @name Setup and finalization of the run */
    //@{
    /** Initialize this object. Called in the run phase just before a run begins. */
    virtual void doinitrun();

    /** Finalize this object. Called in the run phase just after a run has ended. */
    virtual void dofinish();
    //@}
    
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
    tcEventBasePtr eh;
    int mode;

    /* Data handle for ROOT */
    PrtclEvent *pEvent;
    
    /* ROOT, saving data */
    TTree * tree;
    TFile * file;
    
    Timer *timer;
    int timerStep = 1000;
    //@}
private:
    /* The assignment operator is private and must never be called nor implemented. */
    HerwigppTree & operator=(const HerwigppTree &);

}; // Class HerwigppTree


/* Initialization, closing and administrative stuff: */

void HerwigppTree::doinitrun()
{
    AnalysisHandler::doinitrun();
    string fileName = "particles_herwig";
    mode = 0;
    
    /* In a general multithread-case, generate a thread-unique root file name */
    fileName += "_";
    fileName += generator()->runName();
    fileName += ".root";
    
    try {
        size_t pos = fileName.find("jet_");
        string modeName = fileName.substr(17,pos-17);
        if (modeName=="generic") {
            mode = 0;
        } else if (modeName=="di") {
            mode = 1;
        } else if (modeName=="gamma") {
            mode = 2;
        } else if (modeName=="Z") {
            mode = 3;
        }
    } catch (int e) {
        cout << "Invalid mode name: " << e << endl;
    }
    
    /* Setup a root file */
    file = new TFile (fileName.c_str(),"RECREATE");
    if (!file) {
        cout << "Output file could not be created." << endl;
        return;
    }
    file->SetCompressionLevel(1); 
    
    /* Setup a root tree */
    tree = new TTree ("HerwigTree","Herwig++ particle data.");
    if (!tree) {
        cout << "A tree could not be created." << endl;
        return;
    }
    tree->SetAutoSave(100000000);  /* 0.1 GBytes */
    tree->SetCacheSize(10000000);  /* 10 MBytes */
    TTree::SetBranchStyle(1); /* new style */
    
    /* Connect an event handle with the tree */
    pEvent = new PrtclEvent;
    TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
    branch->SetAutoDelete(kFALSE);
    tree->BranchRef();
    
    /* Timing functions */
    timer = new Timer();
    timer->setParams(generator()->N(),timerStep);
    timer->startTiming();
}

void HerwigppTree::dofinish() 
{
    AnalysisHandler::dofinish();
    
    tree->GetCurrentFile();
    tree->AutoSave("Overwrite");
    file->Close();
    
    delete timer;
    cout << "A tree has been written into a .root file" << endl;  
}

/* *** Attention *** This class-description must be correct. */
DescribeNoPIOClass<HerwigppTree,AnalysisHandler>
  describejetanalysisStoreParticles("jetanalysis::HerwigppTree", "../lib/libHerwigppTree.so");

} // namespace jetanalysis

#endif /* STOREPARTICLES_H */
