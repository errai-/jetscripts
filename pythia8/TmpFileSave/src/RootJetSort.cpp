#include "../include/RootJetSort.h"

Int_t RootJetSort::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t RootJetSort::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void RootJetSort::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("px", &px, &b_event_px);
  fChain->SetBranchAddress("py", &py, &b_event_py);
  fChain->SetBranchAddress("pz", &pz, &b_event_pz);
  fChain->SetBranchAddress("e", &e, &b_event_e);
  fChain->SetBranchAddress("id", &id, &b_event_id);
  fChain->SetBranchAddress("particles", &particles, &b_event_particles);
}

void RootJetSort::InitCI(){
  // Create file on which histogram(s) can be saved.
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
}

void RootJetSort::InitFP(){
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
}

void RootJetSort::Insert(double px, double py, double pz, double e, int _status, int _id){
  fastjet::PseudoJet tmpJet(px,py,pz,e);
  Momentum.push_back(tmpJet);
  status.push_back(_status);
  id_store.push_back(_id);
}

void RootJetSort::Clear(){
  Momentum.clear();
  status.clear();
  id_store.clear();
}

void RootJetSort::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void RootJetSort::EventLoop() {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
//    JetStore storaging;
  cout << nentries << endl;
  timer.set_params(nentries,100);

  timer.start_timing();  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
         
    if (jentry!=0&&jentry%100==0){
      timer.print_time();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( !ParticlesToJetsorterInput() ) continue;

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
    inclusiveJets = clustSeq.inclusive_jets( pTMin );    

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    sortedJets    = sorted_by_pt(inclusiveJets);

    // Inspect the amount of jets
    jetMultipl->Fill( sortedJets.size() );

    JetLoop();
    
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
  }
}
    

bool RootJetSort::ParticlesToJetsorterInput(){
  fjInputs.clear();
  
  for (size_t i = 0; i != particles; ++i) {
    fastjet::PseudoJet particleTemp(px[i],py[i],pz[i],e[i]);
    particleTemp.set_user_index( i ); // To access the info of this particle within this event
    fjInputs.push_back( particleTemp );
  }
  
  if (fjInputs.size() == 0) {
    cout << "Error: event with no final state particles" << endl;
    return false;
  }
  return true;
}

void RootJetSort::JetLoop(){
  int counter = 0;
  for (unsigned int i = 0; i < sortedJets.size(); i++) {
    // only count jets that have |eta| < etamax
    if (fabs(sortedJets[i].pseudorapidity()) > etaMax) continue;
    // Inspect only the two leading jets
    if ( counter++ == 2 ) continue;
    // Check that this is a realistic jet
    jetParts = sortedJets[i].constituents();
    if ( jetParts.size() == 1 ) continue;

    // Fill Pt-profile
    ptProfile->Fill( sortedJets[i].pt() );
    // Particle identification
    // Determine whether a jet is dominated by quarks or by gluons
    
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

//     ParticleLoop();
    
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
    if (etSum){ HistFill(i); }
  }
}


void RootJetSort::ParticleLoop(){
  
    piPlus = 0; piMinus = 0;  pi0Gamma = 0; gamma = 0; 
    kaPlus = 0; kaMinus = 0; kSZero = 0; kLZero = 0; 
    proton = 0; aproton = 0; neutron = 0; aneutron = 0;
    lambda0 = 0; sigma = 0; elecmuon = 0;
    others = 0; etSum = 0;
  
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
}

void RootJetSort::HistFill(int i){
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


void RootJetSort::WriteResults(){

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
  TH1D *gq = gluonQuark->ProjectionX("gluonvsquark","");
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
}