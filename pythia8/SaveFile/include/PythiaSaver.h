
#ifndef PYTHIASAVER_H
#define PYTHIASAVER_H

#include <ctime>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

using namespace Pythia8;

using std::cout;
using std::endl;
using std::vector;

// scripts
#include "../include/MinimalEvent.h"
#include "../include/help_functions.h"
#include "../../JetSorter/jetsorter_auxiliary.h"

class PythiaSaver{
private:
  
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  
  // Reweighting for event generation
  PtHatReweightUserHook ptGenReweight;
  int weightedPt = 1; 

  // An object to store the events
  MinimalEvent tmpEvent;
  // A timer to count ETA's
  Timer timer;
  // File to write into
  std::ofstream output;

  int nEvent; // number of events

public:

  PythiaSaver(){}
  
  PythiaSaver(int _nEvent, char *_fileName);
    
  ~PythiaSaver(){
    output.close();
  }
  
  void EventLoop();
  
  void ParticleLoop();
  
};

#endif // PYTHIASAVER_H