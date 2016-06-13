//////////////////////////////////////////////////////////
// A generic main function for fastjet jet clustering   //
// Hannu Siikonen 15.1.2016                             //
//////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>
// Nice libraries from C
#include <cmath>
#include <cassert>

// ROOT, Trees
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

bool main_handle(int argc, char* argv[], string& fullPath, string& treePath,
                 int& mode, int& definition, string& output, string& output2) 
{
    if (argc != 4) {
        cout << "Usage: ./jetanalysis.exe [Standard form input file name] [path - e.g. './'] [Flavour def.]" << endl;
        cout << "Flavour options:" << endl;
        cout << "1: Physics definition" << endl;
        cout << "2: Physics definition with ghost partons" << endl;
        cout << "3: Physics definition with momentum correction" << endl;
        cout << "4: Combination of 2 and 3" << endl;
        cout << "5: Combination of 2 and 10" << endl;
        cout << "6: Historic parton definition" << endl;
        cout << "7: Historic hadron definition" << endl;
        cout << "8: Hadronic definition" << endl;
        cout << "9: Algorithmic definition" << endl;
        cout << "10: Algotirhmic definition with ghost partons" << endl;
        return false;
    }
    
    definition = std::stoi(argv[3]);
    int generator = -1;
    bool beginning = false;
    
    if ( definition<1||definition>10 ) throw std::runtime_error("Flavor options 1-10");
    
    string input = argv[1], tmpStr = "";
    fullPath = argv[2];
    fullPath += input;        
    for (auto i : input) {
        if (i=='_') {
            if (!beginning && tmpStr=="particles") {
                beginning = true;
            } else if (generator == -1) {
                if (tmpStr=="pythia8") {
                    generator = 1;
                } else if (tmpStr=="herwig") {
                    generator = 2;
                } else if (tmpStr=="pythia6") {
                    generator = 3;
                }
                output += i;
                output2 += i;
            } else {
                if (tmpStr=="generic") {
                    mode = 0;
                } else if (tmpStr=="dijet") {
                    mode = 1;
                } else if (tmpStr=="gammajet") {
                    mode = 2;
                } else if (tmpStr=="Zjet") {
                    mode = 3;
                } else if (tmpStr=="ttbarjet") {
                    mode = 4;
                }
                output += i;
                output2 += i;
                if (definition==1) {
                    output += "physics_";
                    output2 += "physics_";
                } else if (definition==2) {
                    output += "ghostphys_";
                    output2 += "ghostphys_";
                } else if (definition==3) {
                    output += "finalphys_";
                    output2 += "finalphys_";
                } else if (definition==4) {
                    output += "finghostphys_";
                    output2 += "finghostphys_";
                } else if (definition==5) {
                    output += "cmsswphys_";
                    output2 += "cmsswphys_";
                } else if (definition==6) {
                    output += "partonhist_";
                    output2 += "partonhist_";
                } else if (definition==7) {
                    output += "historic_";
                    output2 += "historic_";
                } else if (definition==8) {
                    output += "hadronic_";
                    output2 += "hadronic_";
                } else if (definition==9) {
                    output += "algorithmic_";
                    output2 += "algorithmic_";
                } else if (definition==10) {
                    output += "ghostalgo_";
                    output2 += "ghostalgo_";
                }
            }
            tmpStr = "";
        } else {
            tmpStr += i;
            if (beginning) {
                output += i;
                output2 += i;
            }
        }
    }
    if (!beginning || generator==-1 || mode==-1) throw std::runtime_error("Input file uses wrong format");
    
    /* Generator */
    if (generator==1) {
        treePath = "Pythia8Tree";
    } else if (generator == 2) {
        treePath = "HerwigTree";
    } else if (generator == 3) {
        treePath = "Pythia6Tree";
    }
    return true;
}