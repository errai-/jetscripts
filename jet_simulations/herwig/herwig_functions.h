/////////////////////////////////////////////////////////////////////////
// This holds the auxiliary functions and classes that require herwig  //
// Hannu Siikonen 27.3.2015                                            //
/////////////////////////////////////////////////////////////////////////

#ifndef HERWIG_FUNCTIONS
#define HERWIG_FUNCTIONS 

/* Stdlib header file for input and output. */
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "../generic/help_functions.h"

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"

using std::cout;
using std::endl;
using namespace ThePEG;

/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
static int gammaChecker( tPPtr photon ) {
   tParticleVector parents = photon->parents();
   for (tParticleVector::const_iterator pi0t = parents.begin(); 
                              pi0t != parents.end(); ++pi0t) {
      if ( (*pi0t)->id() == 111 ) {
         Lorentz5Momentum momentum = (*pi0t)->momentum();
         
         ParticleVector children = (*pi0t)->children();
         for (ParticleVector::const_iterator g = children.begin(); 
                              g != children.end(); ++g ) {
            momentum -= (*g)->momentum();
         }
         /* Check momentum difference */
         if ( abs(momentum.t()) < 0.001 ) return 1;
      }
   }
   return 0;
}

/* See HadronAndPartonSelector.cc in cmssw, indicates whether a hadron (used for 
 * flavour inspection) is in an excited state or not. This basically checks
 * whether a hadron has a daughter of the same flavour. 
 * The parameter quarkId should be a PDG quark flavour. */
static int isExcitedHadronState(tPPtr part, int quarkId) {
   assert( quarkId>=0 && quarkId<=6 );

   ParticleVector children = part->children();
   for (ParticleVector::const_iterator child = children.begin();
        child != children.end(); ++child ) {
      if ( HadrFuncs::statusCheck( quarkId, (*child)->id() ) ) return 1;
   }
   return 0;
}

#endif // HERWIG_FUNCTIONS