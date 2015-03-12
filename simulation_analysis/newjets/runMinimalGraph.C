#include "TROOT.h"
#include "TChain.h"

runMinimalGraph(){
  TChain *forest = new TChain("JetTree");
  forest->AddFile("jet_storage.root");

  gROOT->ProcessLine(".L JetProcessor.C+");
  JetProcessor handle(forest);

  handle.Loop();
}
