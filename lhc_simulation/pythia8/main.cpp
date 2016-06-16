//////////////////////////////////////////////////////////
// A generic main function for pythia8 event generation //
// Hannu Siikonen 21.8.2015                             //
//////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <stdexcept>
#include <memory>
#include "Pythia8Tree.h"

using std::string;
using std::cout;
using std::endl;
using std::stoi;

int main(int argc, char **argv)
{
    try {
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

        if (mode == 0) {
            P8GenericTree generatorHandle(settings, fileName, mode);
            generatorHandle.EventLoop();
        } else if (mode == 1) {
            P8DijetTree generatorHandle(settings, fileName, mode);
            generatorHandle.EventLoop();
        } else if (mode == 2) {
            P8GammajetTree generatorHandle(settings, fileName, mode);
            generatorHandle.EventLoop();
        } else if (mode == 3) {
            P8ZmumujetTree generatorHandle(settings, fileName, mode);
            generatorHandle.EventLoop();
        } else if (mode == 4) {
            P8ttbarjetTree generatorHandle(settings, fileName, mode);
            generatorHandle.EventLoop();
        }

    } catch (std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << endl;
    }

    return 0;
}
