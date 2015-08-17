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
#include <algorithm>

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
    bool Pythia8ParticleLoop(Pythia&, Event&, PrtclEvent*, const int);
    int GammaChecker(Event&, int);
    int IsExcitedHadronState(Event&, int, int);
    void GhostParticleAdd(PrtclEvent*, Event&, size_t);
    bool HardProcess(Event&, PrtclEvent*, vector<std::size_t>&, std::size_t, const int );

    /* Main loop for storing events
    *  mode:
    *   0 - generic case
    *   1 - standard dijet
    *   2 - gammajet
    *   3 - Zjet
    *   4 - ttbarjet */
    int Pythia8EventLoop(string settings, string fileName, const int mode )
    {
        /* Init pythia by reading settings from a setting file */
        Pythia pythia; Event& event = pythia.event;
        pythia.readFile(settings.c_str()); 
        pythia.init();
        pythia.settings.listChanged();

        int nEvent = pythia.mode("Main:numberOfEvents");

        /* Try to create a file to write */
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

        /* Setup a custom event timer */
        int timerStep = 1000;
        Timer timer; timer.setParams(nEvent,timerStep); timer.startTiming();

        /* Simulation loop: */
        std::size_t ev = 0;
        while (ev != nEvent) {
            if (!pythia.next()) continue;
//             event.list();

            if (Pythia8ParticleLoop(pythia,event,pEvent,mode)) { tree->Fill(); ++ev;}
            pEvent->Clear();

            if (ev%timerStep==0) timer.printTime();
        }
        //pythia.stat();

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

    /* Returns true if event is to be saved */
    bool Pythia8ParticleLoop(Pythia& pythia, Event& event,PrtclEvent* pEvent,const int mode)
    {
        /* Indicators for stable particle indices saved within the hard process */
        vector<std::size_t> skipIndices;
        
        pEvent->fWeight = pythia.info.weight();
        /* Particle loop may save the same particle twice if it is a ghost particle */
        int hardCount = 0;
        for (std::size_t prt = 0; prt!=event.size(); ++prt){

            /* Check for generic ghost particles (uncomment if hadronic definition is used) */
            //GhostParticleAdd(pEvent,event,prt);

            if ( event[prt].statusAbs()==23 ) { 
                if ( event[prt].isParton() ) {
                    /* Gluon or quark from the hard process */
                    ParticleAdd(pEvent,event[prt],3);
                    ++hardCount;
                } else if ( mode==4 && event[prt].idAbs() < 20 ) {
                    /* Leptons, top events */
                    if (!HardProcess( event, pEvent, skipIndices, prt, mode )) return false;
                }
            }
            
            if ( mode>1 && event[prt].statusAbs()==62 ) {
                bool gammaCase = (mode==2 && event[prt].idAbs()==22);
                bool ZCase = (mode==3 && event[prt].idAbs()==23);
                
                // Exclude leptons and uninteresting cases
                if ( (ZCase && skipIndices.size()<2) || (gammaCase && skipIndices.size()<1) ) {
                    if (!HardProcess( event, pEvent, skipIndices, prt, mode )) return false;
                    ++hardCount;
                }
            }
            
            /* Special final-state particles have already been added */
            if ( mode==3 && std::count(skipIndices.begin(), skipIndices.end(), prt)>0 ) {
                continue;
            }
                
            /* Exclude final-state particles counted in the hard process.
               With a generic event type, save pi0 photons with a status 2. */
            if ( event[prt].isFinal() ) {
                if ( (mode==0) && event[prt].id()==22 &&  GammaChecker(event, prt) ) { 
                    ParticleAdd(pEvent,event[prt],2);
                } else {
                    ParticleAdd(pEvent,event[prt],1);
                }
            }
        }
        /* For ttbar-jets, require that one W has gone to leptons and the other to quarks */
        if (mode==4 && skipIndices.size()!=2 ) {
            return false;
        }
        /* Sanity checks */
        if (   (mode<4 && hardCount !=2)
            || (mode==2 && skipIndices.size()!=1 )
            || (mode==3 && skipIndices.size()!=2 ) 
            || (mode==4 && hardCount != 4) ) {
            cout << "Unexpected hard process structure" << endl;
            event.list();
            return false;
        }
        return true;
    }

    bool HardProcess(Event& event, PrtclEvent* pEvent, vector<std::size_t>& skipIds, 
                     std::size_t prt, const int mode )
    {
        if (mode==2) {
            /* For gammajets, find gamma from the hardest subprocess */
            if (event[prt].isFinal()) {
                skipIds.push_back(prt);
                ParticleAdd(pEvent,event[prt],2);
            }
        } else if (mode==3) {
            /* For Zjets, find the final descendant-muons of the Z0 */
            for ( int daughter : event[prt].daughterList() ) {
                if (event[daughter].idAbs()==13) {
                    skipIds.push_back(daughter);
                }
            }
            
            /* Descend to the final muon forms */
            for ( std::size_t i = 0; i < skipIds.size(); ++i ) {
                while (!event[skipIds[i]].isFinal()) {
                    vector<int> mus = event[skipIds[i]].daughterList();
                    for (int daughter : mus) {
                        if (event[daughter].idAbs()==13) {
                            skipIds[i] = daughter; break;
                        }
                    }
                }
                ParticleAdd(pEvent,event[skipIds[i]],2);
            }
        } else if (mode==4) {
            /* For ttbar-events, seek for the final-state leptons */
            /* For a given neutrino expect neutrino etc. (indicated by type) */
            int type = event[prt].idAbs()%2;
            while (!event[prt].isFinal()) {
                vector<int> leptons = event[prt].daughterList();
                for (int daughter : leptons) {
                    if (event[daughter].idAbs()<20 && event[daughter].idAbs()>10) {
                        if (event[daughter].idAbs()%2==type) {
                            prt = daughter; break;
                        }
                    }
                }
            }
            skipIds.push_back(prt);
            ParticleAdd(pEvent,event[prt],2);
        }
        return true;
    }

/////////////////////////////////////////////////////////
// Not in current production but kept in for reference //
/////////////////////////////////////////////////////////

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
}

#endif // PYTHIA8_FUNCTIONS 
