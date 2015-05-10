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
    string fileName = "pythia8_particles_";
    string settings = "pythia8/";
    int choiceId = 1;
    
    if (argc<2 || argc>3) {
        cout << "Usage: ./pythia8.exe [Number of events] (Settings file index) " << endl;
        return 0;
    } 
    if (argc >= 2) {
        nEvent = stoi(argv[1]);
        assert( nEvent > 0 );
    } 
    if (argc == 3) {
        choiceId = stoi(argv[2]);
    }
    switch (choiceId) {
        case 0:
            settings += "pythiaSetting.cmnd";
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
    fileName += std::to_string(nEvent);
    fileName += ".root";
    cout << "Using settings " << settings << endl;
    cout << "Saving to file " << fileName << endl;
    
    /* Create Pythia instance and set it up to generate hard QCD processes
     * above pTHat = 20 GeV for pp collisions at 14 TeV. */
    return pythia8EventLoop(nEvent, settings, fileName, choiceId);
}
