///////////////////////////////////////////////////////////////////////////
// File for auxiliary functions and classes for pythia8 event generation //
// Hannu Siikonen 03.05.2015                                              //
///////////////////////////////////////////////////////////////////////////

#ifndef PYTHIA8_FUNCTIONS
#define PYTHIA8_FUNCTIONS 

/* Stdlib header file for input and output. */
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>

/* Header file to access Pythia 8 program elements. */
#include "Pythia8/Pythia.h"

/* ROOT */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

/* scripts */
#include "../generic/help_functions.h"
#include "../events/PrtclEvent.h"


using namespace Pythia8;
using std::string;


/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
static int gammaChecker( Event &event, int idx ) 
{
    assert( event.size() > idx );
    int mother = event[idx].mother1();
    if ( event[mother].id() != 111 ) return 0;
    double eDifference = abs( event[mother].e() - 
        event[event[mother].daughter1()].e() -
        event[event[mother].daughter2()].e() );
    
    if ( eDifference < 0.001 ) return 1;
    return 0;
}


/* See: HadronAndPartonSelector.cc in CMSSW. Indicates whether a ghost hadron 
 * is in an excited state or not. Checks whether a hadron has a daughter of 
 * the same flavour. Parameter quarkId is a PDG quark flavour. */
static int isExcitedHadronState(Event& event, int idx, int quarkId) 
{
    assert( event.size() > idx );
    assert( quarkId>=0 && quarkId<=6 );

    int dtr1 = event[idx].daughter1(), dtr2 = event[idx].daughter2();   
    if (dtr2 != 0){
        if (dtr1 < dtr2){
            for (int i = dtr1; i <= dtr2; i++) {
                if ( HadrFuncs::statusCheck(quarkId, event[i].id()) ) return 1;
            }
        } else {
            if ( HadrFuncs::statusCheck(quarkId, event[dtr1].id()) ) return 1;
            if ( HadrFuncs::statusCheck(quarkId, event[dtr2].id()) ) return 1;
        }
    } else if (dtr1 != 0) {
        if ( HadrFuncs::statusCheck(quarkId, event[dtr1].id()) ) return 1;
    }
    return 0;
}


static void particleAdd(PrtclEvent* pEvent, Particle& part, int saveStatus)
{
    pEvent->AddPrtcl(part.px(),part.py(),part.pz(),part.e(),part.id(), 
                     part.charge(),saveStatus);
}

static void ghostParticleAdd(PrtclEvent* pEvent, Event& event, size_t prt)
{
    int id = event[prt].id();
    int status = abs( event[prt].status() );
    int ghostStatus = 0;
    
    /* Interesting ghost partons at status codes 71 and 72, ghost hadrons have an id above 100. */
    if (status==71 || status==72) {
        ghostStatus = 4;
    } else if (abs(id) >= 100) {
        /* A hadron may be in all categories -> no 'else' */
        if (HadrFuncs::hasStrange(id) && !isExcitedHadronState(event,prt,3)) {
            ghostStatus = 5; /* s Hadrons */
        }
        if (HadrFuncs::hasCharm(id) && !isExcitedHadronState(event,prt,4)) {  
            ghostStatus = 6; /* c Hadrons */
        }
        if (HadrFuncs::hasBottom(id) && !isExcitedHadronState(event,prt,5)) {
            ghostStatus = 7; /* b Hadrons */
        } 
    }

    if (ghostStatus) { particleAdd(pEvent,event[prt],ghostStatus); }
}

static bool pythia8ParticleLoop(Pythia&, Event&,PrtclEvent*,const int);

/* Main loop for storing events 
 * mode:
 *  0 - generic dijet
 *  1 - standard dijet
 *  2 - gammajet
 *  3 - Zjet */
static int pythia8EventLoop(int nEvent, string settings, string fileName, 
                            const int mode ) 
{
    /* Init pythia */
    Pythia pythia; Event& event = pythia.event;
    pythia.readFile(settings.c_str()); pythia.init();
    pythia.settings.listChanged();
    
    /* Try to create a new file */
    TFile *outFile = new TFile(fileName.c_str(), "NEW");
    if (!outFile->IsOpen()) return 1;
    outFile->SetCompressionLevel(1); /* File is compressed */
    
    /* Create a tree. Autosave every 100 Mb, cache of 10 Mb */ 
    TTree *tree = new TTree("Pythia8Tree","Tree filled with pythia8 data.");
    tree->SetAutoSave(100000000);
    tree->SetCacheSize(10000000);    
    TTree::SetBranchStyle(1); /* New branch style */

    /* Connect an event to the tree */
    PrtclEvent *pEvent = new PrtclEvent();
    TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
    branch->SetAutoDelete(kFALSE);
    tree->BranchRef();

    int timerStep = 100;
    Timer timer; timer.setParams(nEvent,timerStep); timer.startTiming();
    /* Simulation loop: */
    for (size_t ev = 0; ev != nEvent; ++ev) {
        if (!pythia.next()) continue;
        if (ev!=0&&ev%timerStep==0) timer.printTime();
        
        if (pythia8ParticleLoop(pythia,event,pEvent,mode)) tree->Fill();
        pEvent->Clear();
    }

    outFile = tree->GetCurrentFile();
    tree->AutoSave("Overwrite");

    delete pEvent;  pEvent = 0;
    outFile->Close();
    return 0;
}

/* Returns true if event is to be saved */
bool pythia8ParticleLoop(Pythia& pythia, Event& event,PrtclEvent* pEvent,const int mode)
{
    int gammaIdx = -1;
    
    pEvent->fWeight = pythia.info.weight();
    /* Particle loop may save the same particle twice if it is in multiple categories */
    for (size_t prt = 0; prt!=event.size(); ++prt){
        bool hardSubProc = event[prt].status()==-23;
        bool gammaCase = (mode==2 && abs(event[prt].id())==22);
        bool leptonCase = (mode==3 && abs(event[prt].id())==13);
        
        /* Check for generic ghost particles */
        ghostParticleAdd(pEvent,event,prt);
        
        if (gammaCase && gammaIdx==-1 && hardSubProc) {
            gammaIdx = prt;
            while (!event[gammaIdx].isFinal()) {
                if (event[gammaIdx].daughter1() != event[gammaIdx].daughter2()) {
                    /* Skip pair production */
                    return false;
                } else {
                    gammaIdx = event[gammaIdx].daughter1();
                }
            }
            particleAdd(pEvent,event[gammaIdx],2);
        }
            
        /* Final state particles */
        if ( event[prt].isFinal() && prt!=gammaIdx ) {
            if (leptonCase) {
                int probeId = event[prt].mother1(), prevId = prt;
                while (abs(event[probeId].id())==13) {
                    prevId = probeId;
                    probeId = event[probeId].mother1();
                }
                if (event[probeId].id()==23) {
                    particleAdd(pEvent,event[prevId],2);
                } else {
                    particleAdd(pEvent,event[prt],1);
                }
            } else {
                int saveStatus = 1; 
                if ( (mode==0) && event[prt].id()==22 &&  gammaChecker(event, prt) ) { 
                    saveStatus = 2;
                }
                particleAdd(pEvent,event[prt],saveStatus);
            }
        }
        
        if ( hardSubProc && !leptonCase && !gammaCase ) {
            particleAdd(pEvent,event[prt],3);
        }
    }
    return true;
}


#endif // PYTHIA8_FUNCTIONS 
