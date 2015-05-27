#include "TROOT.h"
#include "TChain.h"
#include "JetProcessor.C"
#include <string>

void runMinimalGraph(std::string fileName){
  TChain *forest = new TChain("JetTree");
  forest->AddFile(fileName.c_str());

  JetProcessor handle(forest);

  handle.Loop();
}
