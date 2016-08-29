// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

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

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
#include "../../JetSorter/jetsorter_auxiliary.h"
#include "../include/help_functions.h"
#include "../include/MinimalEvent.h"

using namespace Pythia8;
using std::cout;
using std::endl;
using std::vector;
using std::string;


bool ParticlesToJetsorterInput(MinimalEvent *eventHandler, vector<fastjet::PseudoJet> *fjInputs){
    for (size_t i = 0; i != eventHandler->particles; ++i) {
      fastjet::PseudoJet particleTemp(eventHandler->px[i],eventHandler->py[i],
      eventHandler->pz[i],eventHandler->e[i]);
      particleTemp.set_user_index( i ); // To access the info of this particle within this event
      fjInputs->push_back( particleTemp );
    }
    
    if (fjInputs->size() == 0) {
      std::cout << "Error: event with no final state particles" << std::endl;
      return false;
    }
    return true;
}

int main(int argc, char* argv[]) {

  // ROOT initialization:
  // Create the ROOT application environment.
  TApplication theApp("jet_sorting", &argc, argv);

  int ptBins = 48.;
  const double ptRange[]=
 //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
    //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
    //3637, 3832, 4037};//
    
  // Create file on which histogram(s) can be saved.
  TFile outFile("sortedjets.root", "RECREATE");
  TH1D* ptProfile = new TH1D("pT bins","Pt bins", ptBins, ptRange);
  TH1D* jetMultipl = new TH1D("Jet multiplicity","Jet multiplicity",50,0.,50.);
  TProfile gluonQuark("gq","gq",ptBins,ptRange);
  vector<TH1D*> chargeIndicator;
  chargeIndicator.push_back(new TH1D("gluonjet amount","",150,0,150) );
  chargeIndicator.push_back(new TH1D("quarkjet amount","",150,0,150) );
  chargeIndicator.push_back(new TH1D("gluonjet charge","",30,-15,15) );
  chargeIndicator.push_back(new TH1D("quarkjet charge","",30,-15,15) );
  chargeIndicator.push_back(new TH1D("gluonjet wcharge","",200,-1,1) );
  chargeIndicator.push_back(new TH1D("quarkjet wcharge","",200,-1,1) );
  chargeIndicator.push_back(new TH1D("gluonjet w2charge","",200,-0.8,0.8) );
  chargeIndicator.push_back(new TH1D("quarkjet w2charge","",200,-0.8,0.8) );
  chargeIndicator.push_back(new TH1D("gluonjet w","",250,0,1) );
  chargeIndicator.push_back(new TH1D("quarkjet w","",250,0,1) );

  // Book histograms.
  vector<TProfile*> fractionProfilesGluon;
  vector<TProfile*> fractionProfilesQuark;
  vector<TProfile*> fractionProfilesLQuark;
  vector<TProfile*> fractionProfilesHQuark;
  vector<TProfile*> fractionProfilesAll;
  for (int idx = 0; idx != 16; ++idx){
    std::stringstream tmpString("");
    tmpString << "g" << idx;
    fractionProfilesGluon.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
    tmpString.str("");
    tmpString << "q" << idx;
    fractionProfilesQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
    tmpString.str("");
    tmpString << "lq" << idx;
    fractionProfilesLQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
    tmpString.str("");
    tmpString << "hq" << idx;
    fractionProfilesHQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
    tmpString.str("");
    tmpString << "a" << idx;
    fractionProfilesAll.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
  }
  // ROOT initialization ^

  // Fastjet initialization:
  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range

  // Fastjet analysis - select algorithm and parameters
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;
  // Fastjet initialization ^
  
  // Event loop initialization:
  MinimalEvent eventHandler;
  std::ifstream input;
  input.open("pythia8data.txt", std::ifstream::in);
  size_t iEvent = 0;
  
  // Begin event loop. Generate event; skip if generation aborted.
  while (input.good()) {
    // Reset Fastjet input
    fjInputs.clear();
    // Read an event
    eventHandler.Read(&input);
    
    if ( !ParticlesToJetsorterInput(&eventHandler,&fjInputs) ) continue;

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets( pTMin );
    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    sortedJets    = sorted_by_pt(inclusiveJets);

    // TODO: Not possible without more involved structure of MinimalEvent
//     // Initiate jets with ghost partons
//     for (unsigned int i = 0; i != inspectIndices.size(); ++i){
//       fastjet::PseudoJet particleTemp = event[ inspectIndices[i] ];
//       particleTemp *= pow( 10, -18 );
//       particleTemp.set_user_index( -inspectIndices[i] );
//       fjInputs.push_back( particleTemp );
//     }
//     vector <fastjet::PseudoJet> unsortedGhosts, sortedGhosts;
//     fastjet::ClusterSequence clustSeqGhosts(fjInputs, jetDef);
//     unsortedGhosts = clustSeqGhosts.inclusive_jets( pTMin );
//     sortedGhosts = sorted_by_pt(unsortedGhosts);

    // Inspect the amount of jets
    jetMultipl->Fill( sortedJets.size() );

    int counter = 0;
    for (unsigned int i = 0; i < sortedJets.size(); i++) {
      // only count jets that have |eta| < etamax
      if (fabs(sortedJets[i].pseudorapidity()) > etaMax) continue;
      // Inspect only the two leading jets
      if ( counter++ == 2 ) continue;
      // Check that this is a realistic jet
      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      if ( jetParts.size() == 1 ) continue;

      // Fill Pt-profile
      ptProfile->Fill( sortedJets[i].pt() );
      // Particle identification
      // Determine whether a jet is dominated by quarks or by gluons
      cout << std::setprecision(25);
      // TODO: cannot be used without a more intricate MinimalEvent file
//       vector<int> partonHadronFlavour(2,0);
//       int quarkJetCharge = 0;
//       for (unsigned int j = 0; j != sortedGhosts.size(); ++j) {
//         double dR = deltaR( sortedJets[i].phi(), sortedGhosts[j].phi(), 
//           sortedJets[i].eta(), sortedGhosts[j].eta());
//         if ( dR < 0.1 ) {
//           vector<fastjet::PseudoJet> ghostParts = sorted_by_pt(sortedGhosts[j].constituents());
//           for ( unsigned int k = 0; k != ghostParts.size(); ++k ){
//             if (ghostParts[k].user_index() > 0 ) continue;
//             int idx = -ghostParts[k].user_index();
//             int id = event[ idx ].id();
//             int isHadron = 0;
// 	    // Non-quark or gluon
//             if (id > 9 && id != 21) {
//                 if ( isBottom(id) ){ id = 5;
//                 } else if ( isCharm(id) ){ id = 4;
//                 } else if ( isStrange(id) ){ id = 3;
//                 } else if ( isDown(id) ){ id = 2;
//                 } else if ( isUp(id) ){ id = 1;
//                 } else { id = 0; }
//                 if (!id) { if (isExcitedState( event, idx, id ) ) id = 0; }
//                 isHadron = 1;
//             }
//             if (!id) continue;
//               
//             if (id == 5) { 
//               partonHadronFlavour[isHadron] = 5;
//               if (!isHadron) quarkJetCharge = ChargeSign(id);
//             } else if (id == 4 && partonHadronFlavour[isHadron] != 5) { 
//               partonHadronFlavour[isHadron] = 4;
//               if (!isHadron) quarkJetCharge = ChargeSign(id);
//             } else if (partonHadronFlavour[isHadron] == 0) { 
//               partonHadronFlavour[isHadron] = id;
//               if (!isHadron) quarkJetCharge = ChargeSign(id);
//             }
//           }
//           break;
//         }
//       }
      
      // Initialize energy counters 
      // pi+, pi-, pi0gamma, gamma, k+, k-, ks0, kl0, p+, p-, n, an, l0, s, e/mu, others, tot 
      double piPlus = 0, piMinus = 0,  pi0Gamma = 0, gamma = 0, 
        kaPlus = 0, kaMinus = 0, kSZero = 0, kLZero = 0, 
        proton = 0, aproton = 0, neutron = 0, aneutron = 0,
        lambda0 = 0, sigma = 0, elecmuon = 0,
        others = 0;
      double etSum = 0;

//       double partSum=0, chargSum=0, chargWSum=0, chargW2Sum=0, w2=0;
      for (unsigned int j = 0; j != jetParts.size(); ++j) { 
//         partSum++;
//         chargSum += event[ jetParts[j].user_index() ].charge();
//         chargWSum += event[ jetParts[j].user_index() ].charge()*jetParts[j].perp()/sortedJets[i].perp();
//         chargW2Sum += event[ jetParts[j].user_index() ].charge()*pow(jetParts[j].perp()/sortedJets[i].perp(),2);
//         w2 += pow(jetParts[j].perp()/sortedJets[i].perp(),2);
        double tmpEt = eventHandler.e[ jetParts[j].user_index() ];
        etSum += tmpEt;
        int id = eventHandler.id[ jetParts[j].user_index() ];
        if ( id == 211 ) { 
          piPlus += tmpEt;
        } else if ( id == -211 ) { 
          piMinus+= tmpEt;
        } else if ( id == 22 ) { 
          gamma += tmpEt;
	} else if ( id == 20 ) { // pi0 gamma
	  pi0Gamma += tmpEt;
        } else if ( id == 321 ) { 
          kaPlus += tmpEt;
        } else if ( id == -321 ) { 
          kaMinus += tmpEt;
        } else if ( abs( id ) == 310 ) { 
          kSZero += tmpEt;
        } else if ( abs( id ) == 130 ) { 
          kLZero += tmpEt;
        } else if ( id == 2212 ) { 
          proton += tmpEt;
        } else if ( id == -2212 ) { 
          aproton += tmpEt;
        } else if ( id == 2112 ) { 
          neutron += tmpEt;
        } else if ( id == -2112 ) { 
          aneutron += tmpEt;
        } else if ( abs( id ) == 3122 ) {
          lambda0 += tmpEt;
        } else if ( abs( id ) == 3112 ||
          abs( id ) == 3222 ) {
          sigma += tmpEt;
        } else if ( abs( id ) == 11 || abs( id ) == 13 ) {
          elecmuon += tmpEt;
        } else { 
	  others += tmpEt;        
        }
      }
      // TODO: cannot be used without more evolved MinimalEvent
//       int isHadron = (partonHadronFlavour[1]) ? 1 : 0; 
//       if ( abs(sortedJets[i].eta()) < 1.3 && ( sortedJets[i].perp() > 80 && 
//         sortedJets[i].perp() < 120 ) ){
//         if ( partonHadronFlavour[0]==21 ){
//           chargeIndicator[0]->Fill(partSum);
//           chargeIndicator[2]->Fill(chargSum);
//           chargeIndicator[4]->Fill(chargWSum);
//           chargeIndicator[6]->Fill(chargW2Sum);
//           chargeIndicator[8]->Fill(w2);
//         } else if ( partonHadronFlavour[isHadron] < 9 ){
//           
//           chargeIndicator[1]->Fill(partSum);
//           chargeIndicator[3]->Fill(isHadron ? chargSum : quarkJetCharge*chargSum);
//           chargeIndicator[5]->Fill(isHadron ? chargWSum : quarkJetCharge*chargWSum);
//           chargeIndicator[7]->Fill(isHadron ? chargW2Sum : quarkJetCharge*chargW2Sum);
//           chargeIndicator[9]->Fill(w2);
//         }
//       }

      // Fill the histograms:
      // Fill only gluons or only heavy/light quarks
      if (etSum){
        histFiller( fractionProfilesAll, sortedJets[i].pt(), etSum, piPlus, 
          piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
          aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//         if (partonHadronFlavour[0]==21){
//           gluonQuark.Fill( sortedJets[i].pt(), 1);
//           histFiller( fractionProfilesGluon, sortedJets[i].pt(), etSum, piPlus, 
//             piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
//             aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//         } else if (partonHadronFlavour[isHadron]==4 || partonHadronFlavour[isHadron]==5) {
//           gluonQuark.Fill( sortedJets[i].pt(), 0);
//           histFiller( fractionProfilesHQuark, sortedJets[i].pt(), etSum, piPlus, 
//             piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
//             aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//           histFiller( fractionProfilesQuark, sortedJets[i].pt(), etSum, piPlus, 
//             piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
//             aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//         } else if (partonHadronFlavour[isHadron]==1 || partonHadronFlavour[isHadron]==2 || partonHadronFlavour[isHadron]==3 ){
//           gluonQuark.Fill( sortedJets[i].pt(), 0);
//           histFiller( fractionProfilesLQuark, sortedJets[i].pt(), etSum, piPlus, 
//             piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
//             aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//           histFiller( fractionProfilesQuark, sortedJets[i].pt(), etSum, piPlus, 
//             piMinus, pi0Gamma, kaPlus, kaMinus, kSZero, kLZero, proton,
//             aproton, neutron, aneutron, gamma, lambda0, sigma, elecmuon, others );
//         }
      }
    }
  }

  ptProfile->Write();
  jetMultipl->Write();
  for (unsigned int i = 0; i != fractionProfilesGluon.size(); ++i){
    fractionProfilesGluon[i]->Write();
    fractionProfilesQuark[i]->Write();
    fractionProfilesLQuark[i]->Write();
    fractionProfilesHQuark[i]->Write();
    fractionProfilesAll[i]->Write();
  }

  for (int i = 0; i != 10; ++i){
    chargeIndicator[i]->Write();
  }
  TH1D *gq = gluonQuark.ProjectionX("gluonvsquark","");
  TCanvas *canv = new TCanvas("c1","c1",600,600);
  canv->cd();
  setTDRStyle();
  canv->UseCurrentStyle();
  canv->SetLogx();
  gq->GetXaxis()->SetNoExponent();
  gq->GetXaxis()->SetMoreLogLabels();
  gq->Draw();
  gPad->WaitPrimitive();
  gq->Write();


  // Done.
  return 0;
}
