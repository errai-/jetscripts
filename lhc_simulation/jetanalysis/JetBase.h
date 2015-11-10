/////////////////////////////////////////////////
// The general class for turning particle data //
// into analyzed jet data.                     //
// Hannu Siikonen 7.11.2015                    //
/////////////////////////////////////////////////

#ifndef JETBASE_H
#define JETBASE_H

// General
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>
#include <stdexcept>

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector2.h>

// FastJet interface
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "Pythia8Plugins/FastJet3.h"

// Header file for the classes stored in the TTree if any.
#include "../events/PrtclEvent.h"
#include "../events/JetEvent.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include "../generic/help_functions.h"

using std::cout;
using std::endl;
using std::vector;
using std::cerr;
using fastjet::PseudoJet;

class JetBase
{
public :
    JetBase() : fInitialized(false) {
        cerr << "The default constructor is not intended to be used" << endl;
    }
    JetBase(TTree *, const char*, const char*, Int_t, Int_t);
    
    /* Initializes the tree that is read */
    virtual void        Init(TTree*);
    
    virtual Int_t       GetEntry(Long64_t);
    virtual Long64_t    LoadTree(Long64_t);
    virtual void        Show(Long64_t = -1);
    
    virtual void        EventLoop();
    /* Calculate variables for the newly clustered jets */
    virtual void        JetLoop(Int_t);
    /* Study particle types in the clustered jets */
    virtual void        ParticleLoop();
    virtual inline void EventProcessing(Long64_t);
    virtual inline void PostProcessing() { return; };
    
    virtual void        ParticlesToJetsorterInput();
    /* Event type specific cuts */
    virtual Bool_t      SelectionParams();
    
    /* The best jet flavour definition in terms of robustness. */
    virtual void        PhysicsFlavor(unsigned);
    /* Particle identification the modern "Hadronic definition".
     * Determine whether a jet is dominated by quarks or by gluons.
     * Looping stops when a corresponding jet is found.
     * Hadron flavour is used as a dominating feature.
     * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
     * for further information. */
    virtual void        HadronicFlavor(unsigned);
    /* Algorithmic flavor tagging is somewhat similar to hadronic tagging.
     * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
     * for further information. */
    virtual void        AlgorithmicFlavor(unsigned);
    /* A modernized variant of the physics definition:
     * use ghost particle clustering for the hard process partons. */
    virtual void        PhysClusterFlavor(unsigned);
    /* The hadronic definition without hadrons (or modernized algorithmic deifinition) */
    virtual void        AlgoClusterFlavor(unsigned);

protected: 
///////////
// Fastjet:
///////////
    
    /* IMPORTANT: jet sorting parameters */
    double fR;
    double fMinPT;
    /* end of jet sorting parameters */
    
    int fJetsPerEvent;   /* How many leading jets are stored in a run */
    
    fastjet::JetDefinition   fJetDef;

    /* PseudoJet storages */
    vector<PseudoJet> fJetInputs; 
    vector<PseudoJet> fAuxInputs;
    vector<PseudoJet> fSortedJets; 
    vector<PseudoJet> fJetParts; 

/////////
// Input:
/////////  
    TTree          *fChain;   //! 
    Int_t           fCurrent; //! current Tree number in a fChain
    
    /* Fixed size dimensions of array or collections stored in the TTree if any.
     * If this is too small, Segfaults may follow. */
    static const Int_t  kMaxfPrtcls = 5000;
    
    //PrtclEvent    *event;
    Double_t        fWeight;
    Int_t           fPrtcls_;
    Double_t        fX[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fY[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fZ[kMaxfPrtcls];   //[fPrtcls_]
    Double_t        fT[kMaxfPrtcls];   //[fPrtcls_]
    
    Int_t           fPDGCode[kMaxfPrtcls];   //[fPrtcls_]
    Int_t           fAnalysisStatus[kMaxfPrtcls];   //[fPrtcls_]

//////////
// Output:
//////////
    TFile          *fOutFile;
    TFile          *fOutFile2;
    string          fOutFileName;
    TTree          *fOutTree;
    TBranch        *fJetBranch;
    JetEvent       *fJetEvent;
    
//////////
// Others:
//////////
    /* Energy counters: */ 
    PseudoJet       fPiPlus; 
    PseudoJet       fPiMinus; 
    PseudoJet       fPi0Gamma; 
    PseudoJet       fGamma; 
    PseudoJet       fKaPlus; 
    PseudoJet       fKaMinus; 
    PseudoJet       fKSZero; 
    PseudoJet       fKLZero; 
    PseudoJet       fProton; 
    PseudoJet       fAproton; 
    PseudoJet       fNeutron; 
    PseudoJet       fAneutron;
    PseudoJet       fLambda0; 
    PseudoJet       fSigma; 
    PseudoJet       fElec; 
    PseudoJet       fMuon; 
    PseudoJet       fOthers; 
    PseudoJet       fEtSum;
    
    PseudoJet       fMET;
    
    Timer           fTimer;
    
    /* The user may choose to make cuts manually in the final analysis by setting this to false */
    bool            fJetCuts;
    bool            fParamCuts;
    bool            fParticleStudy;
    bool            fInitialized;
    double          fFlavour;
    
    Int_t           fMode;       /* Event type */
    Int_t           fDefinition; /* Flavour definition */
    Int_t           fBookedParton; /* Monitor partons that are paired with jets */
    
    Int_t           fGammaId;    /* gamma-jets */
    Int_t           fLeptonId;   /* ttbar */
    vector<Int_t>   fLeptonList; /* Z-jets, ttbar */
    vector<Int_t>   fPartonList; /* Physics definition */
    
    JetVariables    fJetVars;
};

#endif // JETBASE_H