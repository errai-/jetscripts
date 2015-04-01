////////////////////////////////////////////////////////////////
// This class sorts simulated jets with the fastjet algorithm.// 
// Hannu Siikonen 11.03.2015                                  //
////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
// #define NDEBUG /* Uncomment to skip asserts */
#include <cassert>
#include "pythia8_functions.h"

using std::string;

int main(int argc, char **argv)
{
   size_t nEvent = 400;
   if (argc > 1) nEvent = atoi(argv[1]); assert( nEvent > 0 );
   bool addLeptons = false;
   
   /* Create Pythia instance and set it up to generate hard QCD processes
    * above pTHat = 20 GeV for pp collisions at 14 TeV. */
   string settings = "pythia8/pythiaSettings.cmnd";
   if (argc > 3) {
      settings = "pythia8/";
      settings += argv[3];
   }

   string fileName = "pythia8_particles.root";
   if (argc > 2) {
      fileName = rootFileName( argv[2] );
   }
   
   if (nEvent > 0) {
      return pythia8EventLoop(nEvent, addLeptons, settings, fileName);
   } else {
      return 1;
   }
}
