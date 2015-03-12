//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 10 15:13:32 2015 by ROOT version 5.34/24
// from TChain JetTree/
//////////////////////////////////////////////////////////

#ifndef JetProcessor_h
#define JetProcessor_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <Math/GenVector/PxPyPzE4D.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxfJets = 100;

class JetProcessor {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //JetEvent        *event;
   Int_t           fJets_;
   Double_t        fJets_fP4_fCoordinates_fX[kMaxfJets];   //[fJets_]
   Double_t        fJets_fP4_fCoordinates_fY[kMaxfJets];   //[fJets_]
   Double_t        fJets_fP4_fCoordinates_fZ[kMaxfJets];   //[fJets_]
   Double_t        fJets_fP4_fCoordinates_fT[kMaxfJets];   //[fJets_]
   Double_t        fJets_fChf[kMaxfJets];   //[fJets_]
   Double_t        fJets_fNhf[kMaxfJets];   //[fJets_]
   Double_t        fJets_fPhf[kMaxfJets];   //[fJets_]
   Double_t        fJets_fElf[kMaxfJets];   //[fJets_]
   Double_t        fJets_fMuf[kMaxfJets];   //[fJets_]
   Double_t        fJets_fChm[kMaxfJets];   //[fJets_]
   Double_t        fJets_fNhm[kMaxfJets];   //[fJets_]
   Double_t        fJets_fPhm[kMaxfJets];   //[fJets_]
   Double_t        fJets_fElm[kMaxfJets];   //[fJets_]
   Double_t        fJets_fMum[kMaxfJets];   //[fJets_]
   Int_t           fJets_fFlavour[kMaxfJets];   //[fJets_]

   // List of branches
   TBranch        *b_event_fJets_;   //!
   TBranch        *b_fJets_fP4_fCoordinates_fX;   //!
   TBranch        *b_fJets_fP4_fCoordinates_fY;   //!
   TBranch        *b_fJets_fP4_fCoordinates_fZ;   //!
   TBranch        *b_fJets_fP4_fCoordinates_fT;   //!
   TBranch        *b_fJets_fChf;   //!
   TBranch        *b_fJets_fNhf;   //!
   TBranch        *b_fJets_fPhf;   //!
   TBranch        *b_fJets_fElf;   //!
   TBranch        *b_fJets_fMuf;   //!
   TBranch        *b_fJets_fChm;   //!
   TBranch        *b_fJets_fNhm;   //!
   TBranch        *b_fJets_fPhm;   //!
   TBranch        *b_fJets_fElm;   //!
   TBranch        *b_fJets_fMum;   //!
   TBranch        *b_fJets_fFlavour;   //!

   JetProcessor(TTree *tree=0);
   virtual ~JetProcessor();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetProcessor_cxx
JetProcessor::JetProcessor(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("JetTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("JetTree","");
      chain->Add("jet_storage.root/JetTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

JetProcessor::~JetProcessor()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetProcessor::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetProcessor::LoadTree(Long64_t entry)
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

void JetProcessor::Init(TTree *tree)
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

   fChain->SetBranchAddress("fJets", &fJets_, &b_event_fJets_);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fX", fJets_fP4_fCoordinates_fX, &b_fJets_fP4_fCoordinates_fX);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fY", fJets_fP4_fCoordinates_fY, &b_fJets_fP4_fCoordinates_fY);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fJets_fP4_fCoordinates_fZ, &b_fJets_fP4_fCoordinates_fZ);
   fChain->SetBranchAddress("fJets.fP4.fCoordinates.fT", fJets_fP4_fCoordinates_fT, &b_fJets_fP4_fCoordinates_fT);
   fChain->SetBranchAddress("fJets.fChf", fJets_fChf, &b_fJets_fChf);
   fChain->SetBranchAddress("fJets.fNhf", fJets_fNhf, &b_fJets_fNhf);
   fChain->SetBranchAddress("fJets.fPhf", fJets_fPhf, &b_fJets_fPhf);
   fChain->SetBranchAddress("fJets.fElf", fJets_fElf, &b_fJets_fElf);
   fChain->SetBranchAddress("fJets.fMuf", fJets_fMuf, &b_fJets_fMuf);
   fChain->SetBranchAddress("fJets.fChm", fJets_fChm, &b_fJets_fChm);
   fChain->SetBranchAddress("fJets.fNhm", fJets_fNhm, &b_fJets_fNhm);
   fChain->SetBranchAddress("fJets.fPhm", fJets_fPhm, &b_fJets_fPhm);
   fChain->SetBranchAddress("fJets.fElm", fJets_fElm, &b_fJets_fElm);
   fChain->SetBranchAddress("fJets.fMum", fJets_fMum, &b_fJets_fMum);
   fChain->SetBranchAddress("fJets.fFlav", fJets_fFlavour, &b_fJets_fFlavour);
   Notify();
}

Bool_t JetProcessor::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetProcessor::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetProcessor::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetProcessor_cxx
