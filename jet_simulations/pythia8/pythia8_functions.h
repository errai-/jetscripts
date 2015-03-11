////////////////////////////////////////////////////////////////////////
// This holds all the auxiliary functions and classes of pythia8 sims //
// Hannu Siikonen 7.3.2015                                            //
////////////////////////////////////////////////////////////////////////

#ifndef PYTHIA8_FUNCTIONS
#define PYTHIA8_FUNCTIONS 

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

#include <cmath>

#include "include/help_functions.h"

using namespace Pythia8;

// A function that checks whether a photon is originated from a pi0 and that
// the energy of the photon-pair corresponds to the pion. returns 0 if
// the origin is not a pion with good energy and 1 if it is
static int gammaChecker( Event &event, int idx ){
   int mother = event[idx].mother1();
   if ( event[mother].id() != 111 ) return 0;
   double eDifference = abs( event[mother].e() -
      event[event[mother].daughter1()].e() -
      event[event[mother].daughter2()].e() );
   
   if ( eDifference < 0.001 ) return 1;
   return 0;
}


static int isExcitedState( Event &event, int idx, int id ) {
   int d1 = event[idx].daughter1(), d2 = event[idx].daughter2();
   if (d2!=0){
      if (d1 < d2){
         for (int i = d1; i <= d2; i++) {
            if ( statusCheck( id, event[i].id() ) ) return 1;
         }
      } else {
         if ( statusCheck( id, event[d1].id() ) ) return 1;
         if ( statusCheck( id, event[d2].id() ) ) return 1;
      }
   } else if (d1!=0) {
      if ( statusCheck( id, event[d1].id() ) ) return 1;
   }
   return 0;
}


#endif // PYTHIA8_FUNCTIONS 
