//////////////////////////////////////////////
// Simple basis for pythia6 functionalities //
// Hannu Siikonen 21.8.2015                 //
//////////////////////////////////////////////

#include "Pythia6Tree.h"

int main(int argc, char** argv)
{
    int nEvent = 400;
    string fileName = "particles_pythia6_";
    int choiceId = 1;
    string nameId = "1";
    int multiplier = 1;

    if (argc<2 || argc>5) {
        cout << "Usage: ./pythia6.exe [Number of events] (mode) (threads)" << endl;
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
            fileName += "generic_";
            break;
        case 1:
            fileName += "dijet_";
            break;
        case 2:
            fileName += "gammajet_";
            break;
        case 3:
            fileName += "Zjet_";
            break;
        case 4:
            fileName += "ttbarjet_";
            break;
        default:
            cout << "Settings file options:" << endl << "0 - generic" << endl
                 << "1 - dijet" << endl << "2 - gammajet" << endl << "3 - Zjet"
                 << "4 - ttbarjet" << endl;
            return 0; break;
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

    try {
        Pythia6Tree generatorHandle(nEvent, fileName, std::stoi(nameId),choiceId);
        generatorHandle.EventLoop();
    } catch (std::exception&e) {
        cout << "An error occurred: " << e.what() << endl;
    }
    return 0;
}
