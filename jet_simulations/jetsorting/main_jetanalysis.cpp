#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
// Nice libraries from C
#include <cmath>
#include <cassert>

// ROOT, Trees
#include "TTree.h"
#include "TChain.h"

// scripts
#include "JetAnalysis.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main(int argc, char* argv[]) {

   /* Mode of operation has to be chosen, since there are a couple of different
    * simulation types in use. Additionally, input and output file names can
    * be specified. */
   if (argc < 2) {
      cout << "Mode of operation: ./a (mode of operation) (read file) (save file)" << endl;
      cout << "Mode of operation has to be entered, either 'p8' for pythia8 or" 
         "'hpp' for herwig++" << endl;
      return 1;
   }
   string treePath;
   string readName;
   string writeName = "jet_storage.root";
   
   /* Mode of operation */
   vector<string> modes;
   modes.push_back("p8"); modes.push_back("hpp");
   if (std::find(modes.begin(),modes.end(),"p8") != modes.end()) {
      treePath = "Pythia8Tree";
      readName = "pythia8_particles.root";
   } else if (std::find(modes.begin(),modes.end(),"hpp") != modes.end()) {
      treePath = "HerwigTree";
      readName = "herwig_particles.root";
   } else {
      cout << "Enter a proper mode of operation (run ./a for further info)" << endl;
      return 1;
   }
   
   /* File to read */
   if (argc > 2) {
      readName = rootFileName( argv[2] );
   }
   
   /* File to write */
   if (argc > 3) {
      writeName = rootFileName( argv[3] );
   }
   
   /* Try to open tree */
   TChain *forest = new TChain(treePath.c_str());
   forest->AddFile(readName.c_str()); /* Opens tree with the highest key */
   if (!(forest->GetEntries())) {
      cout << "Incorrect file name " << readName << " or tree with zero events" << endl;
      return 1;
   }

   /* Analysis process */
   JetAnalysis treeHandle(forest, writeName.c_str());
   treeHandle.EventLoop();
  
   delete forest;
   
   return 0;
}
