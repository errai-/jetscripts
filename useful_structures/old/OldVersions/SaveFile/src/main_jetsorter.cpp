// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

#include "../include/JetSorter.h"

int main(int argc, char* argv[]) {
  
  TApplication theApp("jet_sorting", &argc, argv);

  JetSorter* sorterHandle = new JetSorter();

  sorterHandle->EventLoop();
  
  sorterHandle->WriteResults();

  // Done.
  return 0;
  
}
