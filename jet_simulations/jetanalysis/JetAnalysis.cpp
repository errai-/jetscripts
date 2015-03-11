#include "jetanalysis/JetAnalysis.h"

/////////
// Setup:
/////////

JetAnalysis::JetAnalysis(TTree *tree ) : fChain(0) 
{
   if (tree == 0) {
      TChain * chain = new TChain("Pythia8Tree","");
      chain->Add("particle_storage.root/Pythia8Tree;1");
      tree = chain;
   }
   Init(tree);
  
   jetDef = new fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
  
   // Create a ROOT Tree and one superbranch
   fOutFile = new TFile("jet_storage.root", "RECREATE");
   fOutFile->SetCompressionLevel(1);
   fOutTree = new TTree("JetTree","Tree with jet data");
  
   fjEvent = new JetEvent;
  
   // autosave when 1 Gbyte written 
   fOutTree->SetAutoSave(1000000000);
   // set a 10 MBytes cache (useless when writing local files) 
   fOutTree->SetCacheSize(10000000);  
   TTree::SetBranchStyle(1);
  
   fJetBranch = fOutTree->Branch("event", &fjEvent, 32000, 4);
   fJetBranch->SetAutoDelete(kFALSE);
   fOutTree->BranchRef();
}

JetAnalysis::~JetAnalysis() 
{
   delete fjEvent; fjEvent = 0;
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   if (!fOutTree) return;
   delete fOutTree->GetCurrentFile();
  
   delete jetDef;
}

void JetAnalysis::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fPrtcls", &fPrtcls_, &b_event_fPrtcls_);
   fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fX", 
      fPrtcls_fP4_fCoordinates_fX, &b_fPrtcls_fP4_fCoordinates_fX);
   fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fY", 
      fPrtcls_fP4_fCoordinates_fY, &b_fPrtcls_fP4_fCoordinates_fY);
   fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fZ", 
      fPrtcls_fP4_fCoordinates_fZ, &b_fPrtcls_fP4_fCoordinates_fZ);
   fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fT", 
      fPrtcls_fP4_fCoordinates_fT, &b_fPrtcls_fP4_fCoordinates_fT);
   fChain->SetBranchAddress("fPrtcls.fPDGCode", fPrtcls_fPDGCode, 
                            &b_fPrtcls_fPDGCode);
   fChain->SetBranchAddress("fPrtcls.fChargeTimes3", fPrtcls_fChargeTimes3, 
                            &b_fPrtcls_fChargeTimes3);
   fChain->SetBranchAddress("fPrtcls.IsPi0Photon", fPrtcls_IsPi0Photon, 
                            &b_fPrtcls_IsPi0Photon);
   fChain->SetBranchAddress("fPrtcls.IsJetFlavor", fPrtcls_IsJetFlavor, 
                            &b_fPrtcls_IsJetFlavor);
   fChain->SetBranchAddress("fPrtcls.IsExcitedState", fPrtcls_IsExcitedState,
                            &b_fPrtcls_IsExcitedState);
}

///////////////////////////////////////////////////
// Generic functions, these should not be modified:
///////////////////////////////////////////////////

Int_t JetAnalysis::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t JetAnalysis::LoadTree(Long64_t entry)
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

void JetAnalysis::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

//////////////////////////////////
// Loop over events and particles:
//////////////////////////////////

void JetAnalysis::EventLoop() {
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
  
   // Create file on which a particle data tree is saved (before sampling to jets)
   mTimer.setParams(nentries,500);
   mTimer.startTiming();  
  
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
      mFlavorIndices.clear();
      if (jentry!=0&&jentry%500==0) mTimer.printTime();
    
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
    
      // Initiate jets with ghost partons
      for (size_t i = 0; i != mFlavorIndices.size(); ++i){
         size_t idx = mFlavorIndices[i];
         fastjet::PseudoJet particleTemp(fPrtcls_fP4_fCoordinates_fX[idx],
            fPrtcls_fP4_fCoordinates_fY[idx],fPrtcls_fP4_fCoordinates_fZ[idx],
            fPrtcls_fP4_fCoordinates_fT[idx]);
         particleTemp *= pow( 10, -18 );
         particleTemp.set_user_index( -idx );
         fjInputs.push_back( particleTemp );
      }
      vector <fastjet::PseudoJet> unsortedGhosts;
      fastjet::ClusterSequence clustSeqGhosts(fjInputs, *jetDef);
      unsortedGhosts = clustSeqGhosts.inclusive_jets( pTMin );
      sortedGhosts = sorted_by_pt(unsortedGhosts);

      // Loop through the generated jets
      JetLoop();
    
      fOutTree->Fill();  //fill the tree
      fjEvent->Clear();
   }
  
   fOutFile = fOutTree->GetCurrentFile(); //just in case we switched to a new file
   fOutTree->AutoSave("Overwrite");
   // We own the event (since we set the branch address explicitly), we need to delete it.
   delete fjEvent;  fjEvent = 0;
   
   fOutFile->Close();
}    

bool JetAnalysis::ParticlesToJetsorterInput(){
   fjInputs.clear();
  
   size_t counter = 0;
   for (size_t i = 0; i != fPrtcls_; ++i) {
      if (fPrtcls_IsJetFlavor[i]) { 
         mFlavorIndices.push_back(i);
      } else {
      if (++counter > kMaxfPrtcls){
         cout << "Error: kMaxfPrtcls needs to have a higher value" << endl;
         break;
      }
      fastjet::PseudoJet particleTemp(fPrtcls_fP4_fCoordinates_fX[i],
         fPrtcls_fP4_fCoordinates_fY[i], fPrtcls_fP4_fCoordinates_fZ[i],
         fPrtcls_fP4_fCoordinates_fT[i]);
      particleTemp.set_user_index( i ); // To access the info of this particle within this event
      fjInputs.push_back( particleTemp );
      }
   }
   if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      return false;
   }
   return true;
}

void JetAnalysis::JetLoop(){
   int counter = 0;
   for (size_t i = 0; i < sortedJets.size(); i++) {
      jetParts = sortedJets[i].constituents();
    
      // Some sanity checks/cuts:
      if (fabs(sortedJets[i].pseudorapidity()) > etaMax) continue;
      if ( counter++ == jetsPerEvent ) continue;
      if ( jetParts.size() < 2 ) continue;

      // Find out whether this is a quark or a gluon jet
      FlavorLoop(i);
    
      // Loop over particles within a jet
      ParticleLoop(i);
    
      TypeSort();
    
      fjEvent->AddJet(sortedJets[i].px(),sortedJets[i].py(),sortedJets[i].pz(),
         sortedJets[i].e(),mChf,mNhf,mPhf,mElf,mMuf,mChm,mNhm,mPhm,mElm,mMum,mFlavour);
   }
}


void JetAnalysis::FlavorLoop(size_t i){
  /* Particle identification.
   * Determine whether a jet is dominated by quarks or by gluons.
   * Looping stops when a corresponding jet is found.
   * Hadron flavour is used as a dominating feature.
   * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
   * for further information
   * TODO: set option to determine flavour only based on partons.
   */
   int hadrFlav = 0, partFlav = 0, hardestPartFlav = 0;
   mQuarkJetCharge = 0;
   bool failure = true;
   for (unsigned int j = 0; j != sortedGhosts.size(); ++j) {
      double dR = sortedJets[i].delta_R(sortedGhosts[j]);
      if ( dR < 0.1 ) { // This value could be even smaller
         vector<fastjet::PseudoJet> ghostParts = 
            sorted_by_pt(sortedGhosts[j].constituents());
         for ( unsigned int k = 0; k != ghostParts.size(); ++k ){
            if (ghostParts[k].user_index() > 0 ) continue;
            int idx = -ghostParts[k].user_index();
            int id = fPrtcls_fPDGCode[idx];
            // Hadrons, set the id's to correspond to the hadron flavour
            if (id > 9 && id != 21) {
            if ( isBottom(id) ){
               hadrFlav = 5;
               break;
            } else if ( isCharm(id) && (hadrFlav != 5) ) { 
               hadrFlav = 4;
            }
            //if (!id) { if (fPrtcls_IsExcitedState[idx]) id = 0; }
            } else {
               if (id == 5) { 
                  partFlav = 5;
                  mQuarkJetCharge = ChargeSign(id);
               } else if (id == 4 && partFlav != 5) { 
                  partFlav = 4;
                  mQuarkJetCharge = ChargeSign(id);
               } else if (hardestPartFlav == 0) { 
                  hardestPartFlav = id;
                  mQuarkJetCharge = ChargeSign(id);
               }
            }
         }
      failure = false;
      break;
    }
  }
  if (failure){
    cout << "Error: no matched ghost jet" << endl;
  }
  if (hadrFlav != 0){
    mFlavour = hadrFlav;
  } else {
    mFlavour = TMath::Abs(hardestPartFlav);
  }
}

void JetAnalysis::ParticleLoop(size_t i){
  
  TLorentzVector zero(0,0,0,0);

  mPiPlus = zero; mPiMinus = zero;  mPi0Gamma = zero; mGamma = zero; 
  mKaPlus = zero; mKaMinus = zero; mKSZero = zero; mKLZero = zero; 
  mProton = zero; mAproton = zero; mNeutron = zero; mAneutron = zero;
  mLambda0 = zero; mSigma = zero; mElec = zero, mMuon = zero;
  mOthers = zero; mEtSum = zero;
  
  mPartSum=0; mChargSum=0; mChargWSum=0; mChargW2Sum=0; mW2=0;
  for (unsigned int j = 0; j != jetParts.size(); ++j) { 
    mPartSum++;
    TLorentzVector tmpP( jetParts[j].px(), jetParts[j].py(), jetParts[j].pz(), 
      jetParts[j].e() );

    mChargSum +=  fPrtcls_fChargeTimes3[ jetParts[j].user_index() ]/3.0;
    mChargWSum += (fPrtcls_fChargeTimes3[ jetParts[j].user_index() ]/3.0)
      *jetParts[j].perp()/sortedJets[i].perp();
    mChargW2Sum += (fPrtcls_fChargeTimes3[ jetParts[j].user_index() ]/3.0)
      *pow(jetParts[j].perp()/sortedJets[i].perp(),2);
    mW2 += pow(jetParts[j].perp()/sortedJets[i].perp(),2);
    mEtSum += tmpP;
    int id = fPrtcls_fPDGCode[ jetParts[j].user_index() ];
    if ( id == 211 ) { 
      mPiPlus += tmpP;
    } else if ( id == -211 ) { 
      mPiMinus+= tmpP;
    } else if ( id == 22 ) { 
      mGamma += tmpP;
    } else if ( id == 20 ) { // pi0 gamma
      mPi0Gamma += tmpP;
    } else if ( id == 321 ) { 
      mKaPlus += tmpP;
    } else if ( id == -321 ) { 
      mKaMinus += tmpP;
    } else if ( abs( id ) == 310 ) { 
      mKSZero += tmpP;
    } else if ( abs( id ) == 130 ) { 
      mKLZero += tmpP;
    } else if ( id == 2212 ) { 
      mProton += tmpP;
    } else if ( id == -2212 ) { 
      mAproton += tmpP;
    } else if ( id == 2112 ) { 
      mNeutron += tmpP;
    } else if ( id == -2112 ) { 
      mAneutron += tmpP;
    } else if ( abs( id ) == 3122 ) {
      mLambda0 += tmpP;
    } else if ( abs( id ) == 3112 ||
      abs( id ) == 3222 ) {
      mSigma += tmpP;
    } else if ( abs( id ) == 11 ) {
      mElec += tmpP;
    } else if ( abs( id ) == 13 ) {
      mMuon += tmpP;
    } else { 
      mOthers += tmpP;        
    }
  }    
}

void JetAnalysis::TypeSort()
{
  TLorentzVector zero(0,0,0,0);
  TLorentzVector tmpLorentz = zero;
  
  tmpLorentz += mPiPlus;
  tmpLorentz += mPiMinus;
  tmpLorentz += mKaPlus;
  tmpLorentz += mKaMinus;
  tmpLorentz += mProton;
  tmpLorentz += mAproton;
  tmpLorentz += mSigma;
  mChf = tmpLorentz.Et()/mEtSum.Et();
  mChm = tmpLorentz.M();
  tmpLorentz = zero;
  
  tmpLorentz += mKSZero;
  tmpLorentz += mKLZero;
  tmpLorentz += mNeutron;
  tmpLorentz += mAneutron;
  tmpLorentz += mLambda0;
  mNhf = tmpLorentz.Et()/mEtSum.Et();
  mNhm = tmpLorentz.M();
  tmpLorentz = zero;
  
  tmpLorentz += mPi0Gamma;
  tmpLorentz += mGamma;
  mPhf = tmpLorentz.Et()/mEtSum.Et();
  mPhm = tmpLorentz.M();
  
  mElf = mElec.Et()/mEtSum.Et();
  mElm = mElec.M();
  
  mMuf = mMuon.Et()/mEtSum.Et();
  mMum = mMuon.M();
}
