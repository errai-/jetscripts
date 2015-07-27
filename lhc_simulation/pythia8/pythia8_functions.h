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
#include <vector>

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

namespace
{
    /* A function that checks whether a photon is originated from a pi0 and that
    * the energy of the photon-pair corresponds to the pion. returns 0 if
    * the origin is not a pion with good energy and 1 if it is */
    int GammaChecker( Event &event, int idx ) 
    {
        assert( event.size() > idx );
        
        /* One mother, which is pi0 */
        vector<int> mothers = event[idx].motherList();
        if ( mothers.size()!=1 || event[mothers[0]].id()!=111 ) return 0;
        
        vector<int> daughters = event[mothers[0]].daughterList();
        if ( daughters.size()!=2 ) return 0;
        
        double eDifference = fabs( event[mothers[0]].e() - 
            event[daughters[0]].e() - event[daughters[1]].e() );
        if ( eDifference > 0.001 ) return 0;
        
        return 1;
    }


    /* See: HadronAndPartonSelector.cc in CMSSW. Indicates whether a ghost hadron 
    * is in an excited state or not. Checks whether a hadron has a daughter of 
    * the same flavour. Parameter quarkId is a PDG quark flavour. */
    int IsExcitedHadronState(Event& event, int idx, int quarkId) 
    {
        assert( event.size() > idx );
        assert( quarkId>=0 && quarkId<=6 );

        vector<int> daughters = event[idx].daughterList();
        for (int& dtr : daughters) {
            if ( HadrFuncs::StatusCheck(quarkId, event[dtr].id()) ) return 1;
        }
        return 0;
    }


    bool Pythia8ParticleLoop(Pythia&, Event&,PrtclEvent*,const int);

    /* Main loop for storing events 
    * mode:
    *  0 - generic case
    *  1 - standard dijet
    *  2 - gammajet
    *  3 - Zjet */
    int Pythia8EventLoop(string settings, string fileName, const int mode )
    {
        /* Init pythia with a custom seed */
        Pythia pythia; Event& event = pythia.event;

        pythia.readFile(settings.c_str()); 
        pythia.init();
        pythia.settings.listChanged();
        int nEvent = pythia.mode("Main:numberOfEvents");
        
        /* Try to create a new file */
        TFile *outFile = new TFile(fileName.c_str(), "RECREATE");
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

        int timerStep = 1000;
        Timer timer; timer.setParams(nEvent,timerStep); timer.startTiming();
        /* Simulation loop: */
        std::size_t ev = 0;
        while (ev != nEvent) {
            if (!pythia.next()) continue; //event.list();

            if (Pythia8ParticleLoop(pythia,event,pEvent,mode)) { tree->Fill(); ++ev;}
            if (ev%timerStep==0) timer.printTime();
            pEvent->Clear();
        }

        pythia.stat();

        outFile = tree->GetCurrentFile();
        tree->AutoSave("Overwrite");

        delete pEvent;  pEvent = 0;
        outFile->Close();
        return 0;
    }

    void ParticleAdd(PrtclEvent* pEvent, Particle& part, int saveStatus)
    {
        pEvent->AddPrtcl(part.px(),part.py(),part.pz(),part.e(),part.id(),saveStatus);
    }

    void GhostParticleAdd(PrtclEvent* pEvent, Event& event, size_t prt)
    {
        int id = event[prt].id();
        int status = abs( event[prt].status() );
        int ghostStatus = 0;
        
        /* Interesting ghost partons at status codes 71 and 72, ghost hadrons have an id above 100. */
        if (status==71 || status==72) {
            ghostStatus = 4;
        } else if ( event[prt].isHadron() ) {
            /* A hadron may be in all categories -> no 'else' */
            if (HadrFuncs::HasStrange(id) && !IsExcitedHadronState(event,prt,3)) {
                ghostStatus = 5; /* s Hadrons */
            }
            if (HadrFuncs::HasCharm(id) && !IsExcitedHadronState(event,prt,4)) {  
                ghostStatus = 6; /* c Hadrons */
            }
            if (HadrFuncs::HasBottom(id) && !IsExcitedHadronState(event,prt,5)) {
                ghostStatus = 7; /* b Hadrons */
            }
        }

        if (ghostStatus) { ParticleAdd(pEvent,event[prt],ghostStatus); }
    }

    /* Returns true if event is to be saved */
    bool Pythia8ParticleLoop(Pythia& pythia, Event& event,PrtclEvent* pEvent,const int mode)
    {
        int gammaIdx = -1, mu1Idx = -1, mu2Idx = -1;
        pEvent->fWeight = pythia.info.weight();
        /* Particle loop may save the same particle twice if it is a ghost particle */
        int hardCount = 0;
        for (size_t prt = 0; prt!=event.size(); ++prt){
            bool hardSubProc = event[prt].status()==-23 || event[prt].status()==-22;
            bool gammaCase = (mode==2 && abs(event[prt].id())==22);
            bool ZCase = (mode==3 && abs(event[prt].id())==23);
            
            /* Check for generic ghost particles (consumes space, switched off) */
            //GhostParticleAdd(pEvent,event,prt);
            
            /* Interesting info from the hardest subprocess */
            if (hardSubProc) {
                ++hardCount;
                if (gammaCase) {
                    /* For gammajets, find gamma from the hardest subprocess */
                    gammaIdx = prt;
                    while (!event[gammaIdx].isFinal()) {
                        /* Skip pair production */
                        if (event[gammaIdx].daughterList().size()>1) return false;
                        gammaIdx = event[gammaIdx].daughter1();
                    }
                    ParticleAdd(pEvent,event[gammaIdx],2);
                } else if (ZCase) {

                    /* Note: the muons are sometimes given as an end product for
                    * the hard subprocess, but not always. Z0 is realiable */
                    int probeId = event[prt].daughter1(), prevId = prt;
                    /* Chase down the last Z0 */
                    while (abs(event[probeId].id())==23) {
                        prevId = probeId;
                        probeId = event[probeId].daughter1();
                    }
                    /* Choose the mu-pair coming from Z0 */
                    for ( int daughter : event[prevId].daughterList() ) {
                        if (abs(event[daughter].id())==13) {
                            if (mu1Idx==-1) { mu1Idx = daughter; }
                            else { mu2Idx = daughter; break; }
                        }
                    }
                    if (mu2Idx==-1) return false;
                    /* Find the descendant-muon of mu1 */
                    while (!event[mu1Idx].isFinal()) {
                        vector<int> mu1D = event[mu1Idx].daughterList();
                        for (int daughterI : mu1D) {
                            if (abs(event[daughterI].id())==13) {
                                mu1Idx = daughterI; break;
                            }
                        }
                    }
                    ParticleAdd(pEvent,event[mu1Idx],2);
                    
                    /* Find the descendant-muon of mu2 */
                    while (!event[mu2Idx].isFinal()) {
                        vector<int> mu2D = event[mu2Idx].daughterList();
                        for (int daughterI : mu2D) {
                            if (abs(event[daughterI].id())==13) {
                                mu2Idx = daughterI; break;
                            }
                        }
                    }
                    ParticleAdd(pEvent,event[mu2Idx],2);
                } else if (event[prt].isParton()) {
                    /* Gluon or quark from the hard process */
                    ParticleAdd(pEvent,event[prt],3);
                }
                
                /* Check for the presence of the expected hard-proc particles */
                if (hardCount==2) {
                    if (mode==2 && gammaIdx==-1) { cout << "Unexpected fail" << endl; return false; }
                    if (mode==3 && mu2Idx==-1) { cout << "Unexpected fail" << endl; return false; }
                }
            }
            /* Special final-state particles have already been added */
            if (prt==gammaIdx || prt==mu1Idx || prt==mu2Idx) {
                continue;
            }
                
            /* Final state particles, gamma/mu from gamma/Z-jets excluded */
            if ( event[prt].isFinal() ) {
                int saveStatus = 1;
                if ( (mode==0) && event[prt].id()==22 &&  GammaChecker(event, prt) ) saveStatus = 2;
                ParticleAdd(pEvent,event[prt],saveStatus);
            }
        }
        if (hardCount<2 || hardCount>4) {
            cout << "Unexpected behaviour for the hardest subprocess" << endl;
        }
        return true;
    }

}

#endif // PYTHIA8_FUNCTIONS 
