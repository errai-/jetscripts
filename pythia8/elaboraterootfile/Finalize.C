#include <iostream>
#include <TROOT.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <math.h>

//#include <Hannouris/QCDPFJet.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void Finalize(Long64_t loopLimit=10000,
  string dataPath="/home/hannu/pythia8/rootexamples/")
{
  cout << "Processing run.C.." << endl;

  // Declaration of leaf types
  Double_t        px;
  Double_t        py;
  Double_t        pz;
  Double_t        e;
  Double_t        eta;
  Double_t        phi;
  Int_t           jetsinevent;
  Double_t        chf;
  Double_t        nhf;
  Double_t        phf;
  Double_t        elf;
  Double_t        muf;

  // List of branches
  TBranch        *b_px;   //!
  TBranch        *b_py;   //!
  TBranch        *b_pz;   //!
  TBranch        *b_e;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_jetsinevent;   //!
  TBranch        *b_chf;  //!
  TBranch        *b_nhf;  //!
  TBranch        *b_phf;  //!
  TBranch        *b_elf;  //!
  TBranch        *b_muf;  //!

  // These are not usually changed
  string treePath="pythia8";
  gROOT->ProcessLine(".L Hannouris/QCDEvent.cc++");
  gROOT->ProcessLine(".L Hannouris/QCDJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDMET.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDPFJet.cc+");
  gROOT->ProcessLine(".L Hannouris/QCDEventHdr.cc++");
  gROOT->ProcessLine(".L Hannouris/LorentzVector.h++");
  gROOT->ProcessLine(".L Auxiliary.C++");

  TChain forest(treePath.c_str());
  Int_t fCurrent; //!current Tree number in a TChain

  dataPath += "pythiajets.root";
  forest.AddFile(dataPath.c_str());

  forest.SetBranchAddress("px", &px, &b_px);
  forest.SetBranchAddress("py", &py, &b_py);
  forest.SetBranchAddress("pz", &pz, &b_pz);
  forest.SetBranchAddress("e", &e, &b_e);
  forest.SetBranchAddress("eta", &eta, &b_eta);
  forest.SetBranchAddress("phi", &phi, &b_phi);

  forest.SetBranchAddress("jetsinevent", &jetsinevent, &b_jetsinevent);

  forest.SetBranchAddress("chf", &chf, &b_chf);
  forest.SetBranchAddress("nhf", &nhf, &b_nhf);
  forest.SetBranchAddress("phf", &phf, &b_phf);
  forest.SetBranchAddress("elf", &elf, &b_elf);
  forest.SetBranchAddress("muf", &muf, &b_muf);

  TFile outFile("fjpythia.root", "RECREATE" );  
  outFile.mkdir("ak5");
  outFile.cd("ak5");
  TTree outTree("ProcessedTree","Finer output for pythia 8 results");

  QCDEvent store;
  outTree.Branch("events", &store );

  Long64_t nEntries = forest.GetEntries();

  int counter = 0;
  QCDEventHdr tmpHdr;
  Auxiliary pfStore;

  for ( Long64_t idx = 0; idx != forest.GetEntries(); ++idx ){
    forest.LoadTree( idx );
    forest.GetEntry( idx );
    if ( counter == 0 ) {
      counter = jetsinevent;
      tmpHdr.setRun(1);
      tmpHdr.setLumi(1);
      tmpHdr.setEvt(idx);
    }

    pfStore.addValue( px, py, pz, e, chf, nhf, phf, elf, muf, 0.0, 0.0, 1.0 );
    
    --counter;
    if ( counter == 0 ){
      store.setEvtHdr( tmpHdr );
      store.setPFJets( pfStore.values() );
      pfStore.clear();
      outTree.Fill();
    }
  }
  outTree.Write();
}

