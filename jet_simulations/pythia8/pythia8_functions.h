/////////////////////////////////////////////////////////////////////////
// This holds the auxiliary functions and classes that require pythia8 //
// Hannu Siikonen 7.3.2015                                             //
/////////////////////////////////////////////////////////////////////////

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
#include "../generic/help_functions.h"

/* Header file to access Pythia 8 program elements. */
#include "Pythia8/Pythia.h"

// ROOT
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
static int gammaChecker( Event &event, int idx ) {
   assert( event.size() > idx );
   int mother = event[idx].mother1();
   if ( event[mother].id() != 111 ) return 0;
   double eDifference = abs( event[mother].e() - 
      event[event[mother].daughter1()].e() -
      event[event[mother].daughter2()].e() );
   
   if ( eDifference < 0.001 ) return 1;
   return 0;
}


/* See HadronAndPartonSelector.cc in cmssw, indicates whether a hadron (used for 
 * flavour inspection) is in an excited state or not. This basically checks
 * whether a hadron has a daughter of the same flavour. 
 * The parameter quarkId should be a PDG quark flavour. */
static int isExcitedHadronState(Event &event, int idx, int quarkId) {
   assert( event.size() > idx );
   assert( quarkId>=0 && quarkId<=6 );

   int dtr1 = event[idx].daughter1(), dtr2 = event[idx].daughter2();   
   if (dtr2 != 0){
      if (dtr1 < dtr2){
         for (int i = dtr1; i <= dtr2; i++) {
            if ( statusCheck(quarkId, event[i].id()) ) return 1;
         }
      } else {
         if ( statusCheck(quarkId, event[dtr1].id()) ) return 1;
         if ( statusCheck(quarkId, event[dtr2].id()) ) return 1;
      }
   } else if (dtr1 != 0) {
      if ( statusCheck(quarkId, event[dtr1].id()) ) return 1;
   }
   return 0;
}


/* Main loop for storing events */
static int pythia8EventLoop(int nEvent, bool addLeptons, string settings, 
   string fileName ) {
   Pythia pythia;
   Event& event = pythia.event;
   pythia.readFile(settings.c_str());
   pythia.init();
   pythia.settings.listChanged();
  
   /* Create file on which a particle data tree is saved 
    * (before sampling to jets) */
   TFile *outFile = new TFile(fileName.c_str(), "RECREATE");
   outFile->SetCompressionLevel(1); /* File is compressed */
   
   TTree *tree = new TTree("Pythia8Tree","Tree filled with pythia8 data.");
   PrtclEvent *pEvent = new PrtclEvent();
   
   /* Autosave when 1 Gbyte written */ 
   tree->SetAutoSave(1000000000);
   /* Set a 10 MBytes cache */
   tree->SetCacheSize(10000000);
   
   TTree::SetBranchStyle(1); /* New branch style */
   TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
   branch->SetAutoDelete(kFALSE);
   tree->BranchRef();

   /* Simulation loop: */
   Timer timer;
   timer.setParams(nEvent,100);
   timer.startTiming();  
   for (size_t ev = 0; ev != nEvent; ++ev) {
      if (!pythia.next()) continue;
      if (ev!=0&&ev%100==0) timer.printTime();
      
      for (size_t prt = 0; prt!=event.size(); ++prt){
         double status = abs( event[prt].status() );
         int tmpId = event[prt].id();
         
         /* Stable particles */
         if (event[prt].isFinal() && event[prt].isVisible()) {  
            int pi0Gamma = 0; 
            if ((tmpId == 22) && gammaChecker( event, prt )) pi0Gamma = 1;
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(),tmpId, event[prt].charge(), pi0Gamma ? 10 : 1);
         } 
         
         /* Ghost partons/particles. There might be a small overlap with stable hadrons. */
         int ghostStatus = 0;
         if (status==71 || status==72) { 
            ghostStatus = 11; /* Partons */
         } else if (abs(tmpId) >= 100) { /* Status codes below this are not hadrons */
            if (hasBottom(tmpId) && !isExcitedHadronState(event,prt,5)) {
               ghostStatus = 12; /* b Hadrons */
            } /* A hadron may be in both categories -> no 'else' */
            if (hasCharm(tmpId) && !isExcitedHadronState(event,prt,4)) {  
               ghostStatus = 13; /* c Hadrons */
            }
         } else if ( addLeptons && (status==1 || status==2) ) {
            ghostStatus = 13; /* Leptons */
         }
         
         if (ghostStatus) {
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(), tmpId, event[prt].charge(), ghostStatus);
         }
      } 
      tree->Fill();
      pEvent->Clear();
   }

   outFile = tree->GetCurrentFile();
   tree->AutoSave("Overwrite");

   delete pEvent;  pEvent = 0;
   outFile->Close();
   return 0;
}

#endif // PYTHIA8_FUNCTIONS 
