// This program creates pythia8 particle data and saves it to a binary file

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
// Nice libraries from C
#include <cmath>

#include "../include/PythiaSaver.h"


int main(int argc, char* argv[]) {


  // Settings
  int  nEvent = 100;

  if (argc > 1){
    nEvent = atoi(argv[1]);
  }

  PythiaSaver eventSaver(nEvent, "pythia8data.txt");

  eventSaver.EventLoop();

  // Done.
  return 0;
}
