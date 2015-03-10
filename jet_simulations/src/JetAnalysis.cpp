#include "JetAnalysis.h"

/////////
// Setup:
/////////

JetAnalysis::JetAnalysis(TTree *tree ) : fChain(0) {
  if (tree == 0) {
    TChain * chain = new TChain("Pythia8Tree","");
    chain->Add("particle_storage.root/Pythia8Tree;1");
    tree = chain;
  }
  Init(tree);
  
  jetDef = new fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
  
  // Create a ROOT Tree and one superbranch
  outFile = new TFile("jet_storage.root", "RECREATE");
  outFile->SetCompressionLevel(1);
  outTree = new TTree("JetTree","Tree with jet data");
  
  jEvent = new JetEvent;
  
  outTree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
  outTree->SetCacheSize(10000000);  // set a 10 MBytes cache (useless when writing local files)
  TTree::SetBranchStyle(1);
  
  jetBranch = outTree->Branch("event", &jEvent, 32000, 4);
  jetBranch->SetAutoDelete(kFALSE);
  outTree->BranchRef();
}

JetAnalysis::~JetAnalysis() {
  delete jEvent; jEvent = 0;
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  if (!outTree) return;
  delete outTree->GetCurrentFile();
  
  delete jetDef;
}

void JetAnalysis::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("fParts", &fParts_, &b_event_fParts_);
  fChain->SetBranchAddress("fParts.fP4.fCoordinates.fX", fParts_fP4_fCoordinates_fX, &b_fParts_fP4_fCoordinates_fX);
  fChain->SetBranchAddress("fParts.fP4.fCoordinates.fY", fParts_fP4_fCoordinates_fY, &b_fParts_fP4_fCoordinates_fY);
  fChain->SetBranchAddress("fParts.fP4.fCoordinates.fZ", fParts_fP4_fCoordinates_fZ, &b_fParts_fP4_fCoordinates_fZ);
  fChain->SetBranchAddress("fParts.fP4.fCoordinates.fT", fParts_fP4_fCoordinates_fT, &b_fParts_fP4_fCoordinates_fT);
  fChain->SetBranchAddress("fParts.fPDGCode", fParts_fPDGCode, &b_fParts_fPDGCode);
  fChain->SetBranchAddress("fParts.fChargeTimes3", fParts_fChargeTimes3, &b_fParts_fChargeTimes3);
  fChain->SetBranchAddress("fParts.IsPi0Photon", fParts_IsPi0Photon, &b_fParts_IsPi0Photon);
  fChain->SetBranchAddress("fParts.IsJetFlavor", fParts_IsJetFlavor, &b_fParts_IsJetFlavor);
  fChain->SetBranchAddress("fParts.IsExcitedState", fParts_IsExcitedState, &b_fParts_IsExcitedState);
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
  
  timer.set_params(nentries,100);
  timer.start_timing();  
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
    FlavorIndices.clear();
    if (jentry!=0&&jentry%100==0) timer.print_time();
    
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
    for (size_t i = 0; i != FlavorIndices.size(); ++i){
      size_t idx = FlavorIndices[i];
      fastjet::PseudoJet particleTemp(fParts_fP4_fCoordinates_fX[idx],
        fParts_fP4_fCoordinates_fY[idx],fParts_fP4_fCoordinates_fZ[idx],
        fParts_fP4_fCoordinates_fT[idx]);
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
    
    outTree->Fill();  //fill the tree
    jEvent->Clear();
  }
  
  outFile = outTree->GetCurrentFile(); //just in case we switched to a new file
  outTree->AutoSave("Overwrite");
  // We own the event (since we set the branch address explicitly), we need to delete it.
  delete jEvent;  jEvent = 0;
   
  outFile->Close();
}    

bool JetAnalysis::ParticlesToJetsorterInput(){
  fjInputs.clear();
  
  size_t counter = 0;
  for (size_t i = 0; i != fParts_; ++i) {
    if (fParts_IsJetFlavor[i]){ 
      FlavorIndices.push_back(i);
    }else{
      if (++counter > kMaxfParts){
        cout << "Error: kMaxfParts needs to have a higher value" << endl;
        break;
      }
      fastjet::PseudoJet particleTemp(fParts_fP4_fCoordinates_fX[i],fParts_fP4_fCoordinates_fY[i],fParts_fP4_fCoordinates_fZ[i],fParts_fP4_fCoordinates_fT[i]);
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
    // only count jets that have |eta| < etamax
    if (fabs(sortedJets[i].pseudorapidity()) > etaMax) continue;
    // Inspect only the 'jetsPerEvent' leading jets
    if ( counter++ == jetsPerEvent ) continue;
    // Check that this is a sensible jet
    jetParts = sortedJets[i].constituents();
    if ( jetParts.size() == 1 ) continue;

    // Find out whether this is a quark or a gluon jet
    FlavorLoop(i);
    
    // Loop over particles within a jet
    ParticleLoop(i);
    
    TypeSort();
    
    jEvent->Build(sortedJets[i].px(),sortedJets[i].py(),sortedJets[i].pz(),
      sortedJets[i].e(),chf,nhf,phf,elf,muf,chm,nhm,phm,elm,mum,flavour);
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
  quarkJetCharge = 0;
  bool failure = true;
  for (unsigned int j = 0; j != sortedGhosts.size(); ++j) {
    double dR = sortedJets[i].delta_R(sortedGhosts[j]);
    if ( dR < 0.1 ) { // This value could be even smaller
      vector<fastjet::PseudoJet> ghostParts = 
        sorted_by_pt(sortedGhosts[j].constituents());
      for ( unsigned int k = 0; k != ghostParts.size(); ++k ){
        if (ghostParts[k].user_index() > 0 ) continue;
        int idx = -ghostParts[k].user_index();
        int id = fParts_fPDGCode[idx];
        // Hadrons, set the id's to correspond to the hadron flavour
        if (id > 9 && id != 21) {
          if ( isBottom(id) ){
            hadrFlav = 5;
            break;
          } else if ( isCharm(id) && (hadrFlav != 5) ) { 
            hadrFlav = 4;
          }
          //if (!id) { if (fParts_IsExcitedState[idx]) id = 0; }
        } else {
          if (id == 5) { 
            partFlav = 5;
            quarkJetCharge = ChargeSign(id);
          } else if (id == 4 && partFlav != 5) { 
            partFlav = 4;
            quarkJetCharge = ChargeSign(id);
          } else if (hardestPartFlav == 0) { 
            hardestPartFlav = id;
            quarkJetCharge = ChargeSign(id);
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
    flavour = hadrFlav;
  } else if ((partFlav == 4) || (partFlav == 5)){
    flavour = hardestPartFlav;
  } else {
    flavour = hardestPartFlav;
  }
}

void JetAnalysis::ParticleLoop(size_t i){
  
  TLorentzVector zero(0,0,0,0);

  piPlus = zero; piMinus = zero;  pi0Gamma = zero; gamma = zero; 
  kaPlus = zero; kaMinus = zero; kSZero = zero; kLZero = zero; 
  proton = zero; aproton = zero; neutron = zero; aneutron = zero;
  lambda0 = zero; sigma = zero; elec = zero, muon = zero;
  others = zero; etSum = zero;
  
  partSum=0; chargSum=0; chargWSum=0; chargW2Sum=0; w2=0;
  for (unsigned int j = 0; j != jetParts.size(); ++j) { 
    partSum++;
    TLorentzVector tmpP( jetParts[j].px(), jetParts[j].py(), jetParts[j].pz(), 
      jetParts[j].e() );

    chargSum +=  fParts_fChargeTimes3[ jetParts[j].user_index() ]/3.0;
    chargWSum += (fParts_fChargeTimes3[ jetParts[j].user_index() ]/3.0)
      *jetParts[j].perp()/sortedJets[i].perp();
    chargW2Sum += (fParts_fChargeTimes3[ jetParts[j].user_index() ]/3.0)
      *pow(jetParts[j].perp()/sortedJets[i].perp(),2);
    w2 += pow(jetParts[j].perp()/sortedJets[i].perp(),2);
    etSum += tmpP;
    int id = fParts_fPDGCode[ jetParts[j].user_index() ];
    if ( id == 211 ) { 
      piPlus += tmpP;
    } else if ( id == -211 ) { 
      piMinus+= tmpP;
    } else if ( id == 22 ) { 
      gamma += tmpP;
    } else if ( id == 20 ) { // pi0 gamma
      pi0Gamma += tmpP;
    } else if ( id == 321 ) { 
      kaPlus += tmpP;
    } else if ( id == -321 ) { 
      kaMinus += tmpP;
    } else if ( abs( id ) == 310 ) { 
      kSZero += tmpP;
    } else if ( abs( id ) == 130 ) { 
      kLZero += tmpP;
    } else if ( id == 2212 ) { 
      proton += tmpP;
    } else if ( id == -2212 ) { 
      aproton += tmpP;
    } else if ( id == 2112 ) { 
      neutron += tmpP;
    } else if ( id == -2112 ) { 
      aneutron += tmpP;
    } else if ( abs( id ) == 3122 ) {
      lambda0 += tmpP;
    } else if ( abs( id ) == 3112 ||
      abs( id ) == 3222 ) {
      sigma += tmpP;
    } else if ( abs( id ) == 11 ) {
      elec += tmpP;
    } else if ( abs( id ) == 13 ) {
      muon += tmpP;
    } else { 
      others += tmpP;        
    }
  }    
}

void JetAnalysis::TypeSort()
{
  TLorentzVector zero(0,0,0,0);
  TLorentzVector tmpLorentz = zero;
  
  tmpLorentz += piPlus;
  tmpLorentz += piMinus;
  tmpLorentz += kaPlus;
  tmpLorentz += kaMinus;
  tmpLorentz += proton;
  tmpLorentz += aproton;
  tmpLorentz += sigma;
  chf = tmpLorentz.Et()/etSum.Et();
  chm = tmpLorentz.M();
  tmpLorentz = zero;
  
  tmpLorentz += kSZero;
  tmpLorentz += kLZero;
  tmpLorentz += neutron;
  tmpLorentz += aneutron;
  tmpLorentz += lambda0;
  nhf = tmpLorentz.Et()/etSum.Et();
  nhm = tmpLorentz.M();
  tmpLorentz = zero;
  
  tmpLorentz += pi0Gamma;
  tmpLorentz += gamma;
  phf = tmpLorentz.Et()/etSum.Et();
  phm = tmpLorentz.M();
  
  elf = elec.Et()/etSum.Et();
  elm = elec.M();
  
  muf = muon.Et()/etSum.Et();
  mum = muon.M();
}
