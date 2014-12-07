// simulates pT-observed/pT_actual in a detector

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

// FastJet interface
#include "Pythia8/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

#include "TRandom.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod12.C"

using namespace Pythia8;

using std::vector;
using std::cout;
using std::endl;

// From CMSSW
class PtHatReweightUserHook : public UserHooks
{
  public:
  PtHatReweightUserHook(double _pt = 15, double _power = 4.5) :
  pt(_pt), power(_power) {}
  virtual ~PtHatReweightUserHook() {}
 
  virtual bool canBiasSelection() { return true; }
 
  virtual double biasSelectionBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent)
  {
    //the variable selBias of the base class should be used;
    if ((sigmaProcessPtr->nFinal() == 2)) {
    selBias = pow(phaseSpacePtr->pTHat() / pt, power);
    return selBias;
    }
    selBias = 1.;
    return selBias;
  }
 
  private:
  double pt, power;
};

// A function that checks whether a photon is originated from a pi0 and that
// the energy of the photon-pair corresponds to the pion. returns 0 if
// the origin is not a pion with good energy and 1 if it is
int gammaChecker( Event &event, int idx ){
  int mother = event[idx].mother1();
  if ( event[mother].id() != 111 ) return 0;
  double eDifference = abs( event[mother].e() - 
    event[event[mother].daughter1()].e() - event[event[mother].daughter2()].e() );
  if ( eDifference < 0.001 ) return 1;
  return 0;
}

double deltaR( double phi1, double phi2, double eta1, double eta2 ){
  double dPhi = phi1 - phi2;
  double dEta = eta1 - eta2;
  return pow( pow( dPhi, 2 ) + pow( dEta, 2 ) , 0.5 );
}

// Main loop begins

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("event_generation", &argc, argv);
  int weightedPt = 1;
  // Settings
  int  nEvent = 10000;

  if (argc>1) nEvent = atoi(argv[1]);

  // Old parameters for calorimeter energy cuts
  // const double a = 1.236;
  // const double b = -0.915;
  // const double c = -0.2;
  // New parameters for HCAL energy cuts
  //const double p0 = 1.08063;
  //const double p1 = -0.129273;
  //const double p2 = -0.130000;
  //const double p3 = -2.40708;
  //const double p4 = 3.23;

  // Fractions that describe the fractions of nh energy - these may
  // be changed to depend on the neutron energy
  double fHCAL = 0.7;
  double fECAL = 1-fHCAL;
  double C = 1;

  //TF1 *fp1 = new TF1("fp1","([0]+[1]*pow(max(x,[4]),[2]))"
  //  " * (max(x,[4])+[3])/max(x,[4])",0,4000);
  TF1 *fp1 = new TF1("fp1","([0]+[1]*pow(x,[2]))");
  fp1->SetParameters( 0.90, -8.74, -2.14 );
  //fp1->SetParameters( 1.08063, -0.129273, -0.130000,  -2.40708, 3.23);

  int ptBins = 48.;
  const double ptRange[]=
 //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
    //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
    //3637, 3832, 4037};//

  const int neutrals[]={2112,310,130,3122,3212,3322,111}; 
  const int neutrSize = 7; 
  // neutron, kszero, klzero, lambda0, sigma0, xi0, pi0
  const int chargeds[]={211,321,2212,3222,3112};
  const int chargSize = 5;
  // pi+-,KK+-,p+-,s+,s-

  int power     = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  // Reweighting for event generation
  PtHatReweightUserHook ptGenReweight;

  if (weightedPt){
    pythia.setUserHooksPtr( &ptGenReweight );
  }

  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  //pythia.particleData.listAll();

  pythia.init( 2212, 2212, 14000.);
  pythia.settings.listChanged();
  double minFrac = 0.8;
  double maxFrac = 1.1;

  // Create file on which histogram(s) can be saved.
  TFile outFile("ptcut.root", "RECREATE");
  TProfile* ptProfile = new TProfile("pT bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* RHCAL = new TProfile("RHCAL","",ptBins,ptRange);
  TProfile* hcalProfile = new TProfile("hcal bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* gev3Profile = new TProfile("3gev bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* ptCaloProfile = new TProfile("pT calo bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* hcalCaloProfile = new TProfile("hcal calo bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* gev3CaloProfile = new TProfile("3gev calo bins","", ptBins, ptRange, minFrac, maxFrac);
  TProfile* hiEProfile = new TProfile("hie bins","", ptBins, ptRange, minFrac, maxFrac);

  // Fastjet analysis - select algorithm and parameters
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);

  // Fastjet input
  vector <fastjet::PseudoJet> fullInputs;
  vector <fastjet::PseudoJet> neutCutInputs;
  vector <fastjet::PseudoJet> neut3GeVInputs;
  vector <fastjet::PseudoJet> neutHCALInputs;
  vector <fastjet::PseudoJet> caloCutInputs;
  vector <fastjet::PseudoJet> calo3GeVInputs;
  vector <fastjet::PseudoJet> caloHCALInputs;
  vector <fastjet::PseudoJet> hiEInputs;

  TRandom randomize;

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Reset Fastjet input
    fullInputs.resize(0);
    neutCutInputs.resize(0);
    neut3GeVInputs.resize(0);
    neutHCALInputs.resize(0);
    caloCutInputs.resize(0);
    calo3GeVInputs.resize(0);
    caloHCALInputs.resize(0);
    hiEInputs.resize(0);
    // Particle loop
    for (int i = 0; i != event.size(); ++i) if (event[i].isFinal()) {
      if ( !event[i].isVisible() ) continue;
      fastjet::PseudoJet particleTmp = event[i];
      particleTmp.set_user_index( i );
      fastjet::PseudoJet unalteredParticle = event[i];
      unalteredParticle.set_user_index( i );
      fullInputs.push_back( particleTmp );
      // Cuts are applied to simulate a detector
      int absId = abs(event[i].id());
      int isNHCALStuff = std::find(neutrals,neutrals+neutrSize,absId)!=(neutrals+neutrSize);
      int isCHCALStuff = std::find(chargeds,chargeds+chargSize,absId)!=(chargeds+chargSize); 
      if ( isNHCALStuff || isCHCALStuff ){ // Neutral hadrons and charged hadrons
        
        double tmpE = particleTmp.e();
        double tmpP = pow( pow(particleTmp.px(), 2) + pow(particleTmp.py(), 2) +
          pow(particleTmp.pz(), 2), 0.5 );
        double coeff = particleTmp.e()/tmpP;

        double RHCALf = fp1->Eval(max(tmpP,7.0));
        //tmpE*(p0+p1*pow(max(tmpE,p4),p2))*((max(tmpE,p4)+p3)/max(tmpE,p4));
        // double RTOT = ( a + b*pow( tmpE, c ) );
        // double RHCALf = ( a + b*pow( fHCAL*tmpE, c ) )/0.8176; // The old way

        RHCAL->Fill( tmpE*fHCAL, RHCALf );

        fastjet::PseudoJet pseudoPhoton(particleTmp.px()*coeff,particleTmp.py()*coeff,
          particleTmp.pz()*coeff,particleTmp.e());

        if ( isCHCALStuff ){ // Charged hadrons, neutrals are not touched
          neut3GeVInputs.push_back( unalteredParticle );
          neutHCALInputs.push_back( unalteredParticle );
          neutCutInputs.push_back( unalteredParticle );

          if ( tmpE > 100 ){
            // 50 % fail, of which 50 % EH (not good calibr)
            if ( randomize.Uniform() > 0.75 ){
              particleTmp *= 0.7;
            }
          }
        }
        
        hiEInputs.push_back( particleTmp );
        pseudoPhoton *= fECAL;
        particleTmp *= ((tmpE < 3) ? 1 : RHCALf)*fHCAL*C;
        if (isNHCALStuff){
          neutHCALInputs.push_back( particleTmp );
          neutHCALInputs.push_back( pseudoPhoton );
        }  
        caloHCALInputs.push_back( particleTmp );
        caloHCALInputs.push_back( pseudoPhoton );
        if ( tmpE < 3 ) continue;
        if (isNHCALStuff){
          neutCutInputs.push_back( pseudoPhoton );
          neut3GeVInputs.push_back( pseudoPhoton );
        }  
        caloCutInputs.push_back( pseudoPhoton );
        calo3GeVInputs.push_back( pseudoPhoton );
        if (isNHCALStuff){
          neut3GeVInputs.push_back( unalteredParticle*fHCAL );
        }
        calo3GeVInputs.push_back( unalteredParticle*fHCAL );
      } else { 
        neut3GeVInputs.push_back( particleTmp );
        neutHCALInputs.push_back( particleTmp );
        calo3GeVInputs.push_back( particleTmp );
        caloHCALInputs.push_back( particleTmp );
        hiEInputs.push_back( particleTmp );
      }
      if (!isCHCALStuff)
        neutCutInputs.push_back( particleTmp );
      caloCutInputs.push_back( particleTmp ); 
    }
  
    if (neutCutInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveFullJets, sortedFullJets, 
      inclusiveCutJets, sortedCutJets,
      inclusiveHCALJets, sortedHCALJets, inclusive3GeVJets, sorted3GeVJets,
      inclusiveCaloCutJets, sortedCaloCutJets,
      inclusiveCaloHCALJets, sortedCaloHCALJets, inclusiveCalo3GeVJets, sortedCalo3GeVJets,
      inclusiveHiEJets, sortedHiEJets;

    fastjet::ClusterSequence clustSeqFull(fullInputs, jetDef);
    fastjet::ClusterSequence clustSeqCut(neutCutInputs, jetDef);
    fastjet::ClusterSequence clustSeqHCAL(neutHCALInputs, jetDef);
    fastjet::ClusterSequence clustSeq3GeV(neut3GeVInputs, jetDef);
    fastjet::ClusterSequence caloClustSeqCut(caloCutInputs, jetDef);
    fastjet::ClusterSequence caloClustSeqHCAL(caloHCALInputs, jetDef);
    fastjet::ClusterSequence caloClustSeq3GeV(calo3GeVInputs, jetDef);
    fastjet::ClusterSequence clustSeqHiE(hiEInputs, jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveFullJets = clustSeqFull.inclusive_jets( pTMin );
    sortedFullJets    = sorted_by_pt(inclusiveFullJets);

    inclusiveCutJets = clustSeqCut.inclusive_jets( pTMin );
    sortedCutJets    = sorted_by_pt(inclusiveCutJets);

    inclusiveHCALJets = clustSeqHCAL.inclusive_jets( pTMin );
    sortedHCALJets    = sorted_by_pt(inclusiveHCALJets);

    inclusive3GeVJets = clustSeq3GeV.inclusive_jets( pTMin );
    sorted3GeVJets    = sorted_by_pt(inclusive3GeVJets);

    inclusiveCaloCutJets = caloClustSeqCut.inclusive_jets( pTMin );
    sortedCaloCutJets    = sorted_by_pt(inclusiveCaloCutJets);

    inclusiveCaloHCALJets = caloClustSeqHCAL.inclusive_jets( pTMin );
    sortedCaloHCALJets    = sorted_by_pt(inclusiveCaloHCALJets);

    inclusiveCalo3GeVJets = caloClustSeq3GeV.inclusive_jets( pTMin );
    sortedCalo3GeVJets    = sorted_by_pt(inclusiveCalo3GeVJets);

    inclusiveHiEJets = clustSeqHiE.inclusive_jets( pTMin );
    sortedHiEJets    = sorted_by_pt(inclusiveHiEJets);

    for (unsigned int i = 0; i != sortedFullJets.size(); ++i) {
      // General Cut jets
      double sumCutPt = 0, sum3GeVPt = 0, sumHCALPt = 0;
      if ( abs(sortedFullJets[i].eta()) < 1.3 ) continue; // Eta cut, stabilize
      vector<fastjet::PseudoJet> parts = sorted_by_pt(sortedFullJets[i].constituents());
      if (parts.size()<2) {continue;}
      int successCount = 0;
      for ( unsigned int j = 0; j != sortedCutJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedCutJets[j].phi(),
          sortedFullJets[i].eta(), sortedCutJets[j].eta() );
        if ( dR < 0.25 ){
          sumCutPt += sortedCutJets[j].pt();
          //break;
        }
      }
      if ( sumCutPt > 0 )
        ptProfile->Fill( sortedFullJets[i].pt(), sumCutPt/sortedFullJets[i].pt() );
      
      // Jets with energy missing from hcal
      for ( unsigned int j = 0; j != sortedHCALJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedHCALJets[j].phi(),
          sortedFullJets[i].eta(), sortedHCALJets[j].eta() );
        if ( dR < 0.25 ){
          sumHCALPt += sortedHCALJets[j].pt();
        }
      }
      if ( sumHCALPt > 0 )
        hcalProfile->Fill( sortedFullJets[i].pt() , sumHCALPt/sortedFullJets[i].pt() );

      // Jets with a hcal 3 gev cut
      for ( unsigned int j = 0; j != sorted3GeVJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sorted3GeVJets[j].phi(),
          sortedFullJets[i].eta(), sorted3GeVJets[j].eta() );
        if ( dR < 0.25 ){
          sum3GeVPt += sorted3GeVJets[j].pt();
        }
      }
      if (sum3GeVPt > 0)
        gev3Profile->Fill( sortedFullJets[i].pt() , sum3GeVPt/sortedFullJets[i].pt() );
      
      for ( unsigned int j = 0; j != sortedCaloCutJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedCaloCutJets[j].phi(),
          sortedFullJets[i].eta(), sortedCaloCutJets[j].eta() );
        if ( dR < 0.25 ){
          ptCaloProfile->Fill( sortedFullJets[i].pt(), 
            sortedCaloCutJets[j].pt()/sortedFullJets[i].pt() );
          break;
        }
      }
      // Jets with energy missing from hcal, Calo
      for ( unsigned int j = 0; j != sortedCaloHCALJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedCaloHCALJets[j].phi(),
          sortedFullJets[i].eta(), sortedCaloHCALJets[j].eta() );
        if ( dR < 0.25 ){
          hcalCaloProfile->Fill( sortedFullJets[i].pt(), 
            sortedCaloHCALJets[j].pt()/sortedFullJets[i].pt() );
          break;
        }
      }
      // Jets with a hcal 3 gev cut, Calo
      for ( unsigned int j = 0; j != sortedCalo3GeVJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedCalo3GeVJets[j].phi(),
          sortedFullJets[i].eta(), sortedCalo3GeVJets[j].eta() );
        if ( dR < 0.25 ){
          gev3CaloProfile->Fill( sortedFullJets[i].pt(), 
            sortedCalo3GeVJets[j].pt()/sortedFullJets[i].pt() );
          break;
        }
      }
      // Jets with a hiE cut 
      for ( unsigned int j = 0; j != sortedHiEJets.size(); ++j) {
        double dR = deltaR( sortedFullJets[i].phi(), sortedHiEJets[j].phi(),
          sortedFullJets[i].eta(), sortedHiEJets[j].eta() );
        if ( dR < 0.25 ){
          hiEProfile->Fill( sortedFullJets[i].pt(), 
            sortedHiEJets[j].pt()/sortedFullJets[i].pt() );
          break;
        }
      }
    }
  }

  ptProfile->Write();
  hcalProfile->Write();
  gev3Profile->Write();
  ptCaloProfile->Write();
  hcalCaloProfile->Write();
  gev3CaloProfile->Write();
  hiEProfile->Write();

  RHCAL->Write();

  // Done.
  return 0;
}
