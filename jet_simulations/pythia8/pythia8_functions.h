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

/* Header file to access Pythia 8 program elements. */
#include "Pythia8/Pythia.h"
#include "../generic/help_functions.h"

using namespace Pythia8;

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

#endif // PYTHIA8_FUNCTIONS 
