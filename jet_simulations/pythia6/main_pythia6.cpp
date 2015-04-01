//////////////////////////////////////////////
// Simple basis for pythia6 functionalities //
// Hannu Siikonen 1.4.2015                  //
//////////////////////////////////////////////

#include "pythia6_functions.h"


int main(int argc, char** argv)
{
   size_t nEvent = 400;
   if (argc > 1) nEvent = atoi(argv[1]);
   bool addLeptons = false;
   
   if (nEvent > 0) {
      return makeEventSample(nEvent,addLeptons);
   } else {
      return 1;
   }
}
