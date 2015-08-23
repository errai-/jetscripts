//////////////////////////////////////////////////////////
// A generic main function for pythia8 event generation //
// Hannu Siikonen 08.05.2015                            //
//////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <stdexcept>
#include "Pythia8Tree.h"

using std::string;
using std::cout;
using std::endl;
using std::stoi;

int main(int argc, char **argv)
{
    int mode = 1;
    string name = "";
    if (argc<2 || argc>3) {
        cout << "Usage: ./pythia8.exe (mode) (name)" << endl;
        return 0;
    }
    if (argc >= 2) {
        mode = stoi(argv[1]);
    }
    if (argc >= 3) {
        name = argv[2];
    }
    string settings = name;
    settings += ".cmnd";
    string fileName = "particles_pythia8_";
    fileName += name;
    fileName += ".root";

    try {
        Pythia8Tree generatorHandle(settings, fileName, mode);
        generatorHandle.EventLoop();
    } catch (std::exception& e) {
        cout << "An error occurred: " << e.what() << endl;
    }
        
    return 0;
}
