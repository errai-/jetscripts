
#ifndef PTCUT_H
#define PTCUT_H

// General
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <THStack.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TBranch.h>
#include <TLorentzVector.h>

// FastJet interface
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// Header file for the classes stored in the TTree if any.
#include "../events/PrtclEvent.h"
// Fixed size dimensions of array or collections stored in the TTree if any.
#include "../generic/help_functions.h"
// tdrStyle
#include "../generic/tdrstyle_mod1.C"

using std::cout;
using std::endl;
using std::vector;
using std::max;
using std::string;

class PTCut {
private: // Mostly Fastjet and histogramming
  
    Timer timer;
    
    // ROOT initialization:  
    static const int ptBins = 48.;
    const double ptRange[ptBins+1]=
    //{1, 5, 6, 8, 10, 12, 15
        {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
        97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
        507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
        1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
    //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
    //3637, 3832, 4037};//  
  
    // Create the ROOT application environment.
    // Create file on which histogram(s) can be saved.
    TFile* outFile;
    TH1D* jetMultipl;
    TProfile* RHCAL;
    TProfile* ptProfile;
    TProfile* hcalProfile;
    TProfile* gev3Profile;
    TProfile* ptCaloProfile;
    TProfile* hcalCaloProfile;
    TProfile* gev3CaloProfile;
    TProfile* hiEProfile;

    TProfile* histProfile;
    TProfile* histPNh;
    TProfile* histPCh;
    TProfile* histPE;
    TProfile* histPMu;
    TProfile* histPPh;
    
    // Read TProfiles for detector simulations
    TFile* inFile;
    TProfile* pef1;
    TProfile* pfh;
    TProfile* pfe;
    TProfile* pfe0h;
    TProfile* pfe1h;
    
    // Fastjet initialization:
    int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
    double R      = 0.5;    // Jet size.
    double pTMin  = 20.0;   // Min jet pT
    double etaMax = 1.3;    // Pseudorapidity range

    // Fastjet analysis - select algorithm and parameters
    fastjet::JetDefinition* jetDef;

    // Fastjet initialization ^
    
    // Event loop initialization:
    std::ifstream input;
    size_t iEvent = 0;

    vector<fastjet::PseudoJet> sortedJets, jetParts, fjInputs;
    
    // Simulation parameters: amount of energy in HCAL and ECAL
    double fHCAL = 0.7;
    double fECAL = 1-fHCAL;
    
    // Particle codes:
    static const int neutrSize = 7; 
    const int neutrals[neutrSize]={2112,310,130,3122,3212,3322,111}; 
    // neutron, kszero, klzero, lambda0, sigma0, xi0, pi0
    static const int chargSize = 5;
    const int chargeds[chargSize]={211,321,2212,3222,3112};
    // pi+-,KK+-,p+-,s+,s-
  
    // Scaling for the jets
    TF1 *fp1;
    // H hadrons:
    TF1 *fhr;
    // E hadrons:
    TF1 *fer;
    // eH hadrons:
    TF1 *fre;
    // EH hadrons:
    TF1 *frh;

    TLorentzVector sumPPt, sumPHCAL, sumP3GeV, sumPHi, sumCaloPPt, sumCaloPHCAL,
                   sumCaloP3GeV, sumHistPt, histE, histMu, histCh, histNh, histPh;
  
    Int_t operatingMode; // Detector simulation type
  
public : // Interface to the tree and functions
  
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.
    // If this is too small, Segfaults may follow.
    static const Int_t kMaxfParts = 5000;
    // Declaration of leaf types
    //PrtclEvent        *event;
    Int_t           fPrtcls_;
    Double_t        fWeight;
    Double_t        fX[kMaxfParts];   //[fParts_]
    Double_t        fY[kMaxfParts];   //[fParts_]
    Double_t        fZ[kMaxfParts];   //[fParts_]
    Double_t        fT[kMaxfParts];   //[fParts_]
    Int_t           fPDGCode[kMaxfParts];   //[fParts_]
    Int_t           fChargeTimes3[kMaxfParts];   //[fParts_]
    Int_t           fAnalysisStatus[kMaxfParts];   //[fParts_]

    // List of branches
    TBranch        *b_fPrtcls_;   //!
    TBranch        *b_fWeight;   //!
    TBranch        *b_fX;   //!
    TBranch        *b_fY;   //!
    TBranch        *b_fZ;   //!
    TBranch        *b_fT;   //!
    TBranch        *b_fPDGCode;   //!
    TBranch        *b_fChargeTimes3;   //!
    TBranch        *b_fAnalysisStatus;   //!
    
    virtual Int_t    GetEntry(Long64_t);
    virtual void     Init(TTree *tree); // Chain
    virtual void     InitHistos();
    virtual void     InitFuncs();
    virtual void     InitScales();
    virtual void     Show(Long64_t = -1);
    virtual void     EventLoop();
    virtual bool     ParticlesToJetsorterInput();
    virtual void     JetLoop();
    virtual void     ParticleLoop();
    virtual void     HistFill(int);
    virtual void     WriteResults();
    virtual double   HistResp(TLorentzVector, double, double, double, double);
    
    PTCut(TTree *tree, const char* fileName) : fChain(0) 
    {
        Init(tree);
        
        outFile = new TFile(fileName, "RECREATE");
        
        InitFuncs();
        
        InitHistos();
        
        InitScales();
        
        jetDef = new fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
    }

    ~PTCut() {
        if (!fChain) return;
        delete fChain->GetCurrentFile();
        
        delete outFile;
        delete jetDef;
    }
  
};

#endif // PTCUT_H