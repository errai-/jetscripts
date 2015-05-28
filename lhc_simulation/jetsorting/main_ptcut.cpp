#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "PTCut.h"

using std::cout;
using std::endl;
using std::string;


int main(int argc, char* argv[]) {

    if (argc != 2) {
        cout << "Usage: ./jetanalysis.exe [Standard form input file name]" << endl;
        return 1;
    }
    
    int generator = -1;
    bool beginning = false;
    string input = argv[1], tmpStr = "", output = "detector";
    for (auto i : input) {
        if (i=='_') {
            if (!beginning && tmpStr=="particles") {
                beginning = true;
            } else if (generator == -1) {
                if (tmpStr=="pythia8") {
                    generator = 1;
                }
            } else {
                if (tmpStr!="generic") {
                    cout << "Only generic input should be used for detector simulation" << endl;
                    return 0;
                }
            }
            output += i;
            tmpStr = "";
        } else {
            tmpStr += i;
            if (beginning) {
                output += i;
            }
        }
    }
    if (!beginning || generator==-1) {
        cout << "Input file not of the correct type" << endl;
        return 1;
    }
    
    string treePath;
    
    /* Generator */
    if (generator==1) {
        treePath = "Pythia8Tree";
    } else if (generator == 2) {
        treePath = "HerwigTree";
    } else if (generator == 3) {
        treePath = "Pythia6Tree";
    }
    
    /* Try to open tree */
    TChain *forest = new TChain(treePath.c_str());
    forest->AddFile(input.c_str()); /* Opens tree with the highest key */
    if (!(forest->GetEntries())) {
        cout << "Incorrect file name " << input << " or tree with zero events" << endl;
        return 1;
    }

    /* Analysis process */
    PTCut treeHandle(forest, output.c_str());
  
    treeHandle.EventLoop();
    
    treeHandle.WriteResults();
    
    delete forest;
}