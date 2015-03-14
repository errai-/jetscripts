#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
// Nice libraries from C
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>

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

   /* Create a chain that includes the necessary tree
    * There are a couple of different data options available. */
   if (argc < 2) {
      cout << "Mode of operation: ./a (mode of operation) (read file) (save file)" << endl;
      cout << "Mode of operation has to be entered, either 'p8' for pythia8 or" 
         "'hpp' for herwig++" << endl;
      return 1;
   }
   string treePath;
   string fileName;
   
   vector<string> modes;
   modes.push_back("p8"); modes.push_back("hpp");
   if (std::find(modes.begin(),modes.end(),"p8") != modes.end()) {
      treePath = "Pythia8Tree";
      fileName = "pythia8_particles.root";
   } else if (std::find(modes.begin(),modes.end(),"hpp") != modes.end()) {
      treePath = "HerwigTree";
      fileName = "herwig_particles.root";
   } else {
      cout << "Enter a proper mode of operation (run ./a for further info)" << endl;
      return 1;
   }
   
   if (argc > 2) {
      fileName = rootFileName( argv[2] );
   }
   
   cout << argc << endl;
   string outFile = "jet_storage.root";
   if (argc > 3) {
      outFile = rootFileName( argv[3] );
      cout << outFile << endl;
   }
   
   TChain *forest = new TChain(treePath.c_str());
   /* This opens the tree with the highest key value with the given treePath */
   forest->AddFile(fileName.c_str());
 
   if (forest->GetEntries()==0) {
      cout << "Incorrect file name " << fileName << " or tree with zero events" << endl;
      return 1;
   }

   JetAnalysis treeHandle(forest, outFile.c_str());
   treeHandle.EventLoop();
  
   delete forest;
   
   return 0;
}
