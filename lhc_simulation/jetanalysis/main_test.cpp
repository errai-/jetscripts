#include "common_main.h"

int main(int argc, char* argv[]) 
{
    try {
        int mode = -1;
        int definition;
        string fullPath, treePath;
        string output = "", output2 = "hists_";
        if (!main_handle(argc,argv,fullPath,treePath,mode,definition,output,output2))
            return 0;
        TFile* f;
        TTree* tree;
        
        /* Try to open a tree */
        f = TFile::Open(fullPath.c_str(),"READ");
        if (!f || f->IsZombie()) throw std::runtime_error("Error reading file");
        
        tree = (TTree*)f->Get(treePath.c_str());
        if (!tree || tree->IsZombie()) throw std::runtime_error("Tree could not be opened");

        if (!(tree->GetEntries())) throw std::runtime_error("Zero events found");

        /* Analysis process */
        Test treeHandle(tree, output.c_str(), output2.c_str(), mode, definition);
        treeHandle.EventLoop();
        delete tree;
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << endl;
    }
    
    return 0;
}
