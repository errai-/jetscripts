// main92LinkDef.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header used to generate a ROOT dictionary for the PYTHIA classes.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

#ifdef __CINT__
#include <vector>
#include <JetEvent.h>

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class JetEvent+;
#pragma link C++ class JetData+;
#pragma link C++ class std::vector<Int_t>+;
#endif
