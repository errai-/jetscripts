// simulates pT-observed/pT_actual in a detector

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

#include <cmath>
#include <ctime>
// FastJet interface
#include "Pythia8Plugins/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

#include "TRandom.h"
#include "TLorentzVector.h"

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

double HistResp(TLorentzVector tmpP, double respHCAL, double respECAL,
  double respEHHCAL, double respEHECAL,TProfile* pef1,
  TProfile* pfh,TProfile* pfe,TProfile* pfe0h,TProfile* pfe1h){
  double respSum = 0;
  // In case of inconsistensies in binning, the bins are checked explicitely
  int binHadrf = 0;
  int binEcalf = 0;
  double trueFracECAL = 0;
  double fE = 0;
  // H: all energy in HCAL; this is multiplied with HCAL response
  binHadrf = pfh->FindBin(tmpP.Perp());
  respSum += (pfh->GetBinContent(binHadrf))*respHCAL; 
  
  // E: energy only in ECAL, fraction multiplied with ECAL response
  binHadrf = pfe->FindBin(tmpP.Perp());
  respSum += (pfe->GetBinContent(binHadrf))*respECAL; 
  
  // pf0: nothing is added for the lost hadrons, effectively:
  // respSum += (pf0->GetBinContent(binHadrf))*0;
  
  // eH: energy in both ECAL and HCAL, mostly HCAL. Does not work
  // intuitively (compared to EH). First, approximate H by HCAL
  // resp. times tot energy - 0.3 GeV. Move 0.3 GeV to ECAL and
  // use response for this.
  binHadrf = pfe0h->FindBin(tmpP.Perp());
  respSum += (pfe0h->GetBinContent(binHadrf))*( 0.3*(respECAL-1)
    + respHCAL );

  // EH: energy in both ECAL and HCAL
  binHadrf = pfe1h->FindBin(tmpP.Perp());
  binEcalf = pef1->FindBin(tmpP.Perp());

  fE = pef1->GetBinContent(binEcalf); 
  trueFracECAL = fE/((1-fE)*respEHECAL/respEHHCAL+1); 
  respSum += (pfe1h->GetBinContent(binHadrf))*(trueFracECAL*respEHECAL
    +(1-trueFracECAL)*respEHHCAL);

  return respSum;
}

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
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();
  pythia.settings.listChanged();

  // To read:
  TFile inFile("pfhadrons.root", "READ");
  TProfile* pef1 = (TProfile*) inFile.Get("pef1");
  TProfile* pfh = (TProfile*) inFile.Get("pfh");
  TProfile* pfe = (TProfile*) inFile.Get("pfe");
  TProfile* pfe0h = (TProfile*) inFile.Get("pfe0h");
  TProfile* pfe1h = (TProfile*) inFile.Get("pfe1h");

  // To write:
  TProfile* histProfile = new TProfile("hist bins","",ptBins, ptRange);
  TProfile* histPNh = new TProfile("nh bins","",ptBins,ptRange);
  TProfile* histPCh = new TProfile("ch bins","",ptBins,ptRange);
  TProfile* histPE = new TProfile("e bins","",ptBins,ptRange);
  TProfile* histPMu = new TProfile("mu bins","",ptBins,ptRange);
  TProfile* histPPh = new TProfile("ph bins","",ptBins,ptRange);

  TProfile* ptProfile = new TProfile("pT bins","", ptBins, ptRange);
  TProfile* RHCAL = new TProfile("RHCAL","",ptBins,ptRange);
  TProfile* hcalProfile = new TProfile("hcal bins","", ptBins, ptRange);
  TProfile* gev3Profile = new TProfile("3gev bins","", ptBins, ptRange);
  TProfile* ptCaloProfile = new TProfile("pT calo bins","", ptBins, ptRange);
  TProfile* hcalCaloProfile = new TProfile("hcal calo bins","", ptBins, ptRange);
  TProfile* gev3CaloProfile = new TProfile("3gev calo bins","", ptBins, ptRange);
  TProfile* hiEProfile = new TProfile("hie bins","", ptBins, ptRange);
  
  // Response functions:
  
  // Old:
  TF1 *fp1 = new TF1("fp1","([0]+[1]*pow(max(x,[4]),[2]))"
    " * (max(x,[4])+[3])/max(x,[4])",0,4000);
  fp1->SetParameters(1.0806, -0.0838, -0.1300, -2.3051, 3.44);
  // Old forms:
  //fp1->SetParameters( 1.08063, -0.129273, -0.130000,  -2.40708, 3.23);
  //tmpE*(p0+p1*pow(max(tmpE,p4),p2))*((max(tmpE,p4)+p3)/max(tmpE,p4));
  // double RTOT = ( a + b*pow( tmpE, c ) );
  // double RHCALf = ( a + b*pow( fHCAL*tmpE, c ) )/0.8176; // The old way

  // New: // a 3 GeV minimum energy is utilized
  // H hadrons:
  TF1 *fhr = new TF1("fhr","([0]+[1]*pow(max(x,[4])+[3],[2]))" 
    " * (max(x,[4])+[3])/max(x,[4])",0,2000);
  fhr->SetParameters(1.0806,-0.125961,-0.130000,-2.76844,3);
  // E hadrons:
  TF1 *fer = new TF1("fer","([0]+[1]*pow(max(x,[4])+[3],[2]))" 
    " * (max(x,[4])+[3])/max(x,[4])",0,2000);
  fer->SetParameters(1.11506,-2.45363,-0.349480,8.25860,3);
  // eH hadrons:
  TF1 *fre = new TF1("fre","([0]+[1]*pow(max(x,[4])+[3],[2]))" 
    " * (max(x,[4])+[3])/max(x,[4])",0,2000);
  fre->SetParameters(1.35221,-1.13817,-0.1300,-1.79764,3);
  // EH hadrons:
  TF1 *frh = new TF1("frh","([0]+[1]*pow(max(x,[4])+[3],[2]))" 
    " * (max(x,[4])+[3])/max(x,[4])",0,2000);
  frh->SetParameters(1.01771,-0.125961,-0.130000,-3.86871,3);

  // Fastjet analysis - select algorithm and parameters
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);

  // Fastjet input
  vector <fastjet::PseudoJet> fullInputs;

  TRandom randomize;

  std::clock_t start = std::clock();
  double time_processor = 0; int hours; int minutes; int seconds; 
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // event.list();
    if (iEvent!=0&&iEvent%100==0){
      time_processor = (std::clock() - start)/(( (double) CLOCKS_PER_SEC ) );
      time_processor = time_processor*( ((double) nEvent)/iEvent-1);
      minutes =  time_processor/60; hours = minutes/60;
      seconds = time_processor-60*minutes;
      minutes = minutes - hours*60;
      cout << iEvent << " events analyzed. Eta : " << hours << "h" <<
        minutes << "m" << seconds << "s." << endl;
    }
    // Reset Fastjet input
    fullInputs.resize(0);
    // Particle loop
    for (int i = 0; i != event.size(); ++i) if (event[i].isFinal()) {
      if ( !event[i].isVisible() ) continue;
      fastjet::PseudoJet particleTmp = event[i];
      particleTmp.set_user_index( i );
      fullInputs.push_back( particleTmp );
    }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveFullJets, sortedFullJets; 
    fastjet::ClusterSequence clustSeqFull(fullInputs, jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveFullJets = clustSeqFull.inclusive_jets( pTMin );
    sortedFullJets    = sorted_by_pt(inclusiveFullJets);

    for (unsigned int i = 0; i != sortedFullJets.size(); ++i) {
      // General Cut jets
      // Limitations: more than 1 particle in a jet, eta < 1.3
      if (abs(sortedFullJets[i].eta()) > 1.3) continue;
      vector<fastjet::PseudoJet> parts = sorted_by_pt(sortedFullJets[i].constituents());
      if (parts.size()<2) {continue;}
      
      TLorentzVector sumPPt(0,0,0,0), sumPHCAL(0,0,0,0), sumP3GeV(0,0,0,0), sumPHi(0,0,0,0);
      TLorentzVector sumCaloPPt(0,0,0,0), sumCaloPHCAL(0,0,0,0), sumCaloP3GeV(0,0,0,0);
      TLorentzVector sumHistPt(0,0,0,0), histE(0,0,0,0), histMu(0,0,0,0);
      TLorentzVector histCh(0,0,0,0), histNh(0,0,0,0), histPh(0,0,0,0);
      
      for ( unsigned int j = 0; j != parts.size(); ++j ){
        Particle tmpPart = event[parts[j].user_index()]; 
        int absId = abs( tmpPart.id() ); 

        int isNHCALStuff = std::find(neutrals,neutrals+neutrSize,absId)!=(neutrals+neutrSize);
        int isCHCALStuff = std::find(chargeds,chargeds+chargSize,absId)!=(chargeds+chargSize); 
        
        TLorentzVector tmpP( tmpPart.px(), tmpPart.py(), tmpPart.pz(), tmpPart.e() );
        double RHCALf = 1;

        // In case of hadrons, set initial values for temporary containers
        if ( isNHCALStuff || isCHCALStuff ){
          RHCALf = fp1->Eval(max(7.0,tmpP.Perp()));  
          //RHCAL->Fill( tmpP.Perp(), RHCALf ); // Record fractions
        }

        if ( isCHCALStuff ){ // Charged hadrons

          // In 50 % (tracking) of the cases 50 % of charged hadr. receive a worse response (0.7) 
          double cutTrackAndHighPt = 1; 
          if ( tmpP.Perp() > 100 ){
            cutTrackAndHighPt = 0.925;
          } else {
            cutTrackAndHighPt = 0.97;
          }

          double response = (fECAL+fHCAL*((tmpP.Perp()<3) ? 1 : RHCALf) );
          
          double lowECut = ( tmpP.Perp() < 3 ) ? 0 : 1; // 3 GeV cut

          sumHistPt += cutTrackAndHighPt*tmpP;
          histCh += cutTrackAndHighPt*tmpP;

          sumP3GeV += tmpP;
          sumCaloP3GeV += lowECut*tmpP;
          sumPHCAL += tmpP;
          sumCaloPHCAL += response*tmpP;
          sumPPt += cutTrackAndHighPt*tmpP;
          sumCaloPPt += cutTrackAndHighPt*lowECut*response*tmpP;
          sumPHi += cutTrackAndHighPt*tmpP;
        
        } else if ( isNHCALStuff ){ // Neutral hadrons
          double response = (fECAL+fHCAL*((tmpP.Perp()<3) ? 1 : RHCALf) );
          
          double lowECut = ( tmpP.Perp() < 3 ) ? 0 : 1; // 3 GeV cut
          
          double respHCAL = fhr->Eval( tmpP.Perp() );
          double respECAL = fer->Eval( tmpP.Perp() );
          double respEHHCAL = frh->Eval( tmpP.Perp() );
          double respEHECAL = fre->Eval( tmpP.Perp() );
          respHCAL = (respHCAL > 0) ? respHCAL : 0;
          respECAL = (respECAL > 0) ? respECAL : 0;
          respEHHCAL = (respEHHCAL > 0) ? respEHHCAL : 0;
          respEHECAL = (respEHECAL > 0) ? respEHECAL : 0;

          double histResp = HistResp(tmpP, respHCAL, respECAL, respEHHCAL, 
            respEHECAL, pef1, pfh, pfe, pfe0h, pfe1h);
          histResp = ( (histResp > 0) ? histResp : 0 );  
          RHCAL->Fill(tmpP.Perp(),histResp);

          sumHistPt += lowECut*histResp*tmpP;
          histNh += lowECut*histResp*tmpP;

          sumP3GeV += lowECut*tmpP; 
          sumCaloP3GeV += lowECut*tmpP;
          sumPHCAL += response*tmpP;
          sumCaloPHCAL += response*tmpP;
          sumPPt += lowECut*response*tmpP;
          sumCaloPPt += lowECut*response*tmpP;
          sumPHi += tmpP;
        
        } else { // Others

          sumHistPt += tmpP;
          if ( absId == 22 ){
            histPh += tmpP;
          } else if ( absId == 11 ){
            histE += tmpP;
          } else if ( absId == 13 ){
            histMu += tmpP;
          }

          sumP3GeV += tmpP; 
          sumCaloP3GeV += tmpP;
          sumPHCAL += tmpP; 
          sumCaloPHCAL += tmpP;
          sumPPt += tmpP;
          sumCaloPPt += tmpP; 
          sumPHi += tmpP;
        
        }
      }
      histProfile->Fill( sortedFullJets[i].pt(), sumHistPt.Perp()/sortedFullJets[i].pt() );
      histPNh->Fill( sortedFullJets[i].pt(), histNh.Perp()/sumHistPt.Perp() );
      histPCh->Fill( sortedFullJets[i].pt(), histCh.Perp()/sumHistPt.Perp() );
      histPE->Fill( sortedFullJets[i].pt(), histE.Perp()/sumHistPt.Perp() );
      histPMu->Fill( sortedFullJets[i].pt(), histMu.Perp()/sumHistPt.Perp() );
      histPPh->Fill( sortedFullJets[i].pt(), histPh.Perp()/sumHistPt.Perp() );

      ptProfile->Fill( sortedFullJets[i].pt(), sumPPt.Perp()/sortedFullJets[i].pt() );
      hcalProfile->Fill( sortedFullJets[i].pt(), sumPHCAL.Perp()/sortedFullJets[i].pt() );
      gev3Profile->Fill( sortedFullJets[i].pt(), sumP3GeV.Perp()/sortedFullJets[i].pt() );
      ptCaloProfile->Fill( sortedFullJets[i].pt(), sumCaloPPt.Perp()/sortedFullJets[i].pt() );
      hcalCaloProfile->Fill( sortedFullJets[i].pt(), sumCaloPHCAL.Perp()/sortedFullJets[i].pt() );
      gev3CaloProfile->Fill( sortedFullJets[i].pt(), sumCaloP3GeV.Perp()/sortedFullJets[i].pt() );
      hiEProfile->Fill( sortedFullJets[i].pt(), sumPHi.Perp()/sortedFullJets[i].pt() );
    }  
  
    if (fullInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }
  }
  // Create file on which histogram(s) can be saved.
  TFile outFile("ptcut_hists_full.root", "RECREATE");
  
  ptProfile->Write();
  hcalProfile->Write();
  gev3Profile->Write();
  ptCaloProfile->Write();
  hcalCaloProfile->Write();
  gev3CaloProfile->Write();
  hiEProfile->Write();
  cout << "fuufuu\n";
      
  histProfile->Write();
  histPNh->Write();
  histPCh->Write();
  histPE->Write();
  histPMu->Write();
  histPPh->Write();

  RHCAL->Write();
  // Done.
  return 0;
}
