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
    int choiceId = 1;
    string name = "";
    if (argc<2 || argc>3) {
        cout << "Usage: ./pythia8.exe (mode) (name)" << endl;
        return 0;
    }
    if (argc >= 2) {
        choiceId = stoi(argv[1]);
    }
    if (argc >= 3) {
        name = argv[2];
    }
    string settings = name;
    settings += ".cmnd";
    string fileName = "particles_pythia8_";
    fileName += name;
    fileName += ".root";
    
    return Pythia8EventLoop(settings, fileName, choiceId);
}
