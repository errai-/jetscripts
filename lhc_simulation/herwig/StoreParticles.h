///////////////////////////////////////////////////////////////////
// A class for Herwig to store simulation data into a .root file //
// Based on a template: contains some peculiarities              //
// Hannu Siikonen 5.6.2015                                       //
///////////////////////////////////////////////////////////////////

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

class StoreParticles: public AnalysisHandler {

public:

    StoreParticles() {}
    
    ~StoreParticles() {}
    
    /** Analyze a given Event.
        * @param event pointer to the Event to be analyzed.
        * @param ieve the event number.
        * @param loop the number of times this event has been presented.
        * @param state nonzero if the event has been manipulated. */
    virtual void analyze(tEventPtr event, long ieve, int loop, int state);


    /** Standard Init function, called exactly once. */
    static void Init() { static ClassDocumentation<StoreParticles> documentation("No documentation"); }

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
    TTree * herwigTree;
    TFile * herwigFile;
    //@}
private:
    /* The assignment operator is private and must never be called nor implemented. */
    StoreParticles & operator=(const StoreParticles &);
};

}

#endif /* STOREPARTICLES_H */
