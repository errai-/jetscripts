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


int main(int argc, char* argv[]) 
{
    if (argc != 3) {
        cout << "Usage: ./jetanalysis [Standard form input file name] [Flavour def.]" << endl;
        cout << "Flavour options:" << endl << "1: Physics definition" << endl;
        cout << "2: Hadronic definition" << endl;
        return 1;
    }
    
    int definition = std::stoi(argv[2]);
    if ( definition!=1&&definition!=2 ) {
        cout << "Flavour options: 1/2" << endl;
        return 0;
    }
    
    int generator = -1, mode = -1;
    bool beginning = false;
    string input = argv[1], tmpStr = "", output = "jets_";
    for (auto i : input) {
        if (i=='_') {
            if (!beginning && tmpStr=="particles") {
                beginning = true;
            } else if (generator == -1) {
                if (tmpStr=="pythia8") {
                    generator = 1;
                }
                output += i;
            } else {
                if (tmpStr=="generic") {
                    mode = 0;
                } else if (tmpStr=="dijet") {
                    mode = 1;
                } else if (tmpStr=="gammajet") {
                    mode = 2;
                } else if (tmpStr=="Zjet") {
                    mode = 3;
                }
                output += i;
            }
            tmpStr = "";
        } else {
            tmpStr += i;
            if (beginning) {
                output += i;
            }
        }
    }
    if (!beginning || generator==-1 || mode==-1) {
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
    JetAnalysis treeHandle(forest, output.c_str(), mode, definition);
    treeHandle.EventLoop();
    
    delete forest;
    
    return 0;
}
