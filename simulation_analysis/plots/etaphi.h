//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 13 16:56:17 2016 by ROOT version 6.05/01
// from TTree JetTree/Tree with jet data
// found on file: pythia8_dijet_historic_100.root
//////////////////////////////////////////////////////////

#ifndef EtaPhi_h
#define EtaPhi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/PxPyPzE4D.h"

class EtaPhi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxfJets = 5000;

   // Declaration of leaf types
 //JetEvent        *event;
   Int_t           fJets;
   Double_t        fX[kMaxfJets];   //[fJets_]
   Double_t        fY[kMaxfJets];   //[fJets_]
   Double_t        fZ[kMaxfJets];   //[fJets_]
   Double_t        fT[kMaxfJets];   //[fJets_]
   Int_t           fFlav[kMaxfJets];   //[fJets_]
   Float_t         fChf[kMaxfJets];   //[fJets_]
   Float_t         fNhf[kMaxfJets];   //[fJets_]
   Float_t         fPhf[kMaxfJets];   //[fJets_]
   Float_t         fElf[kMaxfJets];   //[fJets_]
   Float_t         fMuf[kMaxfJets];   //[fJets_]
   Float_t         fChm[kMaxfJets];   //[fJets_]
   Float_t         fNhm[kMaxfJets];   //[fJets_]
   Float_t         fPhm[kMaxfJets];   //[fJets_]
   Float_t         fElm[kMaxfJets];   //[fJets_]
   Float_t         fMum[kMaxfJets];   //[fJets_]
   Float_t         fPartonPT[kMaxfJets];   //[fJets_]
   Float_t         fMatchPT[kMaxfJets];   //[fJets_]
   Int_t           fConstituents[kMaxfJets];   //[fJets_]
   Float_t         fPTD[kMaxfJets];   //[fJets_]
   Float_t         fSigma2[kMaxfJets];   //[fJets_]
   Float_t         fDR[kMaxfJets];   //[fJets_]
   Float_t         fDR_Next[kMaxfJets];   //[fJets_]
   Float_t         fAlpha[kMaxfJets];   //[fJets_]
   Float_t         fDPhi[kMaxfJets];   //[fJets_]
   Double_t        fWeight;

   EtaPhi(TTree *tree=0);
   virtual ~EtaPhi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EtaPhi_cxx
EtaPhi::EtaPhi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("target.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("target.root");
      }
      f->GetObject("JetTree",tree);

   }
   Init(tree);
}

EtaPhi::~EtaPhi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EtaPhi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EtaPhi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EtaPhi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fJets", &fJets);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
   fChain->SetBranchAddress("fJets.fFlav", fFlav);
   fChain->SetBranchAddress("fJets.fChf", fChf);
   fChain->SetBranchAddress("fJets.fNhf", fNhf);
   fChain->SetBranchAddress("fJets.fPhf", fPhf);
   fChain->SetBranchAddress("fJets.fElf", fElf);
   fChain->SetBranchAddress("fJets.fMuf", fMuf);
   fChain->SetBranchAddress("fJets.fChm", fChm);
   fChain->SetBranchAddress("fJets.fNhm", fNhm);
   fChain->SetBranchAddress("fJets.fPhm", fPhm);
   fChain->SetBranchAddress("fJets.fElm", fElm);
   fChain->SetBranchAddress("fJets.fMum", fMum);
   fChain->SetBranchAddress("fJets.fPartonPT", fPartonPT);
   fChain->SetBranchAddress("fJets.fMatchPT", fMatchPT);
   fChain->SetBranchAddress("fJets.fConstituents", fConstituents);
   fChain->SetBranchAddress("fJets.fPTD", fPTD);
   fChain->SetBranchAddress("fJets.fSigma2", fSigma2);
   fChain->SetBranchAddress("fJets.fDR", fDR);
   fChain->SetBranchAddress("fJets.fDR_Next", fDR_Next);
   fChain->SetBranchAddress("fJets.fAlpha", fAlpha);
   fChain->SetBranchAddress("fJets.fDPhi", fDPhi);
   fChain->SetBranchAddress("fWeight", &fWeight);
   Notify();
}

Bool_t EtaPhi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EtaPhi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EtaPhi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EtaPhi_cxx
