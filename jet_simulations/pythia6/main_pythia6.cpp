
#include "pythia6_functions.h"


int main(int argc, char** argv)
{

  Int_t n = 100;
  if (argc > 1)
    n = strtol(argv[1], NULL, 0);

  int retVal = 0;
  if (n > 0)
    retVal = makeEventSample(n);
  else {
    retVal = showEventSample();
  }

  return retVal;
}
