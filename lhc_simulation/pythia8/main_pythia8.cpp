//////////////////////////////////////////////////////////
// A generic main function for pythia8 event generation //
// Hannu Siikonen 08.05.2015                            //
//////////////////////////////////////////////////////////

#include <iostream>
#include <string>
// #define NDEBUG /* Uncomment to skip asserts */
#include <cassert>
#include "pythia8_functions.h"

using std::string;
using std::cout;
using std::endl;
using std::stoi;

int main(int argc, char **argv)
{
    int nEvent = 400;
    string fileName = "particles_pythia8_";
    string settings = "";
    int choiceId = 1;
    string nameId = "1";
    int multiplier = 1;
    
    if (argc<2 || argc>5) {
        cout << "Usage: ./pythia8.exe [Number of events] (Settings file index) (Number of threads to be used)" << endl;
        return 0;
    } 
    if (argc >= 2) {
        nEvent = stoi(argv[1]);
        assert( nEvent > 0 );
    } 
    if (argc >= 3) {
        choiceId = stoi(argv[2]);
    }
    if (argc >= 4) {
        nameId = argv[3];
    }
    if (argc >= 5) {
        multiplier = stoi(argv[4]);
    }
    switch (choiceId) {
        case 0:
            settings += "pythia_generic.cmnd";
            fileName += "generic_";
            break;
        case 1:
            settings += "pythia_dijet.cmnd";
            fileName += "dijet_";
            break;
        case 2:
            settings += "pythia_gammajet.cmnd";
            fileName += "gammajet_";
            break;
        case 3:
            settings += "pythia_Zjet.cmnd";
            fileName += "Zjet_";
            break;
        default:
            cout << "Settings file options:" << endl << "0 - generic" 
                    << endl << "1 - dijet" << endl << "2 - gammajet" 
                    << endl << "3 - Zjet" << endl;
            return 0;
            break;
    }
    string fileNameFinal = fileName;
    fileName += std::to_string(nEvent);
    fileNameFinal += std::to_string(nEvent*multiplier);
    fileNameFinal += ".root";
    
    if (multiplier > 1) {
        fileName += "_";
        fileName += nameId;
        if (nameId=="1") {
            TFile *outFile = new TFile(fileNameFinal.c_str(), "RECREATE");
            outFile->Close();
        }
    }
    fileName += ".root";
    
    cout << "Using settings " << settings << endl;
    cout << "Saving to file " << fileName << endl;
    
    return Pythia8EventLoop(nEvent, settings, fileName, choiceId, stoi(nameId));
}
