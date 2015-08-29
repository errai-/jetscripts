//////////////////////////////////////////////////////////
// A generic main function for fastjet jet clustering   //
// Hannu Siikonen 21.8.2015                             //
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

// scripts
#include "JetAnalysis.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main(int argc, char* argv[]) 
{
    if (argc != 4) {
        cout << "Usage: ./jetanalysis.exe [Standard form input file name] [path - e.g. './'] [Flavour def.]" << endl;
        cout << "Flavour options:" << endl << "1: Physics definition" << endl;
        cout << "2: Hadronic definition" << endl;
        return 0;
    }
    
    int definition = std::stoi(argv[3]);
    int generator = -1, mode = -1;
    bool beginning = false;
    TFile *f;
    TTree *tree;
    string output = "", output2 = "hists_";   
    try {
    
        if ( definition!=1&&definition!=2 ) throw std::runtime_error("Flavor options 1/2");
        
        string input = argv[1], fullPath = argv[2], tmpStr = "";
        fullPath += input;        
        string treePath;
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
                    } else {
                        output += "hadronic_";
                        output2 += "hadronic_";
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
        
        /* Try to open a tree */
        f = TFile::Open(fullPath.c_str(),"READ");
        if (!f || f->IsZombie()) throw std::runtime_error("Error reading file");
        
        tree = (TTree*)f->Get(treePath.c_str());
        if (!tree || tree->IsZombie()) throw std::runtime_error("Tree could not be opened");
            
        if (!(tree->GetEntries())) throw std::runtime_error("Zero events found");
        
    } catch (std::exception& e) {
        cout << "An error occurred: " << e.what() << endl;
    }

    /* Analysis process */
    try {
        JetAnalysis treeHandle(tree, output.c_str(), output2.c_str(), mode, definition);
        treeHandle.EventLoop();
        delete tree;
    } catch (std::exception& e) {
        cout << "An error occurred: " << e.what() << endl;
    }
    
    return 0;
}
