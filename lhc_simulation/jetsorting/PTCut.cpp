#include "PTCut.h"

Int_t PTCut::GetEntry(Long64_t entry)
{
// Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

void PTCut::Init(TTree *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("fWeight", &fWeight, &b_fWeight);
    fChain->SetBranchAddress("fPrtcls", &fPrtcls_, &b_fPrtcls_);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fX", fX, &b_fX);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fY", fY, &b_fY);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fZ", fZ, &b_fZ);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fT", fT, &b_fT);
    fChain->SetBranchAddress("fPrtcls.fPDGCode", fPDGCode, &b_fPDGCode);
    fChain->SetBranchAddress("fPrtcls.fChargeTimes3", fChargeTimes3, &b_fChargeTimes3);
    fChain->SetBranchAddress("fPrtcls.fAnalysisStatus", fAnalysisStatus, &b_fAnalysisStatus);    
}

void PTCut::InitHistos()
{
    jetMultipl = new TH1D("Jet multiplicity","Jet multiplicity",50,0.,50.);
        
    RHCAL = new TProfile("RHCAL","",ptBins,ptRange);
    ptProfile = new TProfile("pT bins","", ptBins, ptRange);
    hcalProfile = new TProfile("hcal bins","", ptBins, ptRange);
    gev3Profile = new TProfile("3gev bins","", ptBins, ptRange);
    ptCaloProfile = new TProfile("pT calo bins","", ptBins, ptRange);
    hcalCaloProfile = new TProfile("hcal calo bins","", ptBins, ptRange);
    gev3CaloProfile = new TProfile("3gev calo bins","", ptBins, ptRange);
    hiEProfile = new TProfile("hie bins","", ptBins, ptRange);
    
    histProfile = new TProfile("hist bins","",ptBins, ptRange);
    histPNh = new TProfile("nh bins","",ptBins,ptRange);
    histPCh = new TProfile("ch bins","",ptBins,ptRange);
    histPE = new TProfile("e bins","",ptBins,ptRange);
    histPMu = new TProfile("mu bins","",ptBins,ptRange);
    histPPh = new TProfile("ph bins","",ptBins,ptRange); 
}

void PTCut::InitFuncs()
{
    // a 3 GeV minimum energy is utilized
    fp1 = new TF1("fp1","([0]+[1]*pow(TMath::Max(x,[4]),[2]))"
        " * (TMath::Max(x,[4])+[3])/TMath::Max(x,[4])",0,4000);
    fp1->SetParameters(1.0806, -0.0838, -0.1300, -2.3051, 3.44);
    // Old forms:
    // fp1->SetParameters( 1.08063, -0.129273, -0.130000,  -2.40708, 3.23);
    // tmpE*(p0+p1*pow(max(tmpE,p4),p2))*((max(tmpE,p4)+p3)/max(tmpE,p4));
    // double RTOT = ( a + b*pow( tmpE, c ) );
    // double RHCALf = ( a + b*pow( fHCAL*tmpE, c ) )/0.8176;
    
    fhr  = new TF1("fhr","([0]+[1]*pow(TMath::Max(x,[4])+[3],[2]))" 
        " * (TMath::Max(x,[4])+[3])/TMath::Max(x,[4])",0,2000);
    fhr->SetParameters(1.0806,-0.125961,-0.130000,-2.76844,3);
    
    fer = new TF1("fer","([0]+[1]*pow(TMath::Max(x,[4])+[3],[2]))" 
        " * (TMath::Max(x,[4])+[3])/TMath::Max(x,[4])",0,2000);  
    fer->SetParameters(1.11506,-2.45363,-0.349480,8.25860,3);
    
    fre = new TF1("fre","([0]+[1]*pow(TMath::Max(x,[4])+[3],[2]))" 
        " * (TMath::Max(x,[4])+[3])/TMath::Max(x,[4])",0,2000);
    fre->SetParameters(1.35221,-1.13817,-0.1300,-1.79764,3);
    
    frh = new TF1("frh","([0]+[1]*pow(TMath::Max(x,[4])+[3],[2]))" 
        " * (TMath::Max(x,[4])+[3])/TMath::Max(x,[4])",0,2000); 
    frh->SetParameters(1.01771,-0.125961,-0.130000,-3.86871,3);
}

void PTCut::InitScales()
{
    inFile = new TFile("generic/pfhadrons.root", "READ");
    pef1 = (TProfile*) inFile->Get("pef1");
    pfh = (TProfile*) inFile->Get("pfh");
    pfe = (TProfile*) inFile->Get("pfe");
    pfe0h = (TProfile*) inFile->Get("pfe0h");
    pfe1h = (TProfile*) inFile->Get("pfe1h"); 
}

void PTCut::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

void PTCut::EventLoop() 
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    cout << nentries << endl;
    timer.setParams(nentries,500);

    timer.startTiming();  
    for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
            
        if (jentry!=0&&jentry%500==0) timer.printTime();
        
        fChain->GetEntry(jentry);
        if ( !ParticlesToJetsorterInput() ) continue;

        // Run Fastjet algorithm
        vector <fastjet::PseudoJet> inclusiveJets;
        fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
        inclusiveJets = clustSeq.inclusive_jets( pTMin );    

        // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
        sortedJets    = sorted_by_pt(inclusiveJets);

        // Inspect the amount of jets
        jetMultipl->Fill( sortedJets.size() );

        // Loop through the generated jets
        JetLoop();
    }
}
    

bool PTCut::ParticlesToJetsorterInput()
{
    fjInputs.clear();
    
    for (size_t i = 0; i != fPrtcls_; ++i) {
        if (fAnalysisStatus[i]!=1 && fAnalysisStatus[i]!=2) continue;
        fastjet::PseudoJet particleTemp(fX[i],fY[i],fZ[i],fT[i]);
        particleTemp.set_user_index( abs(fPDGCode[i]) ); 
        fjInputs.push_back( particleTemp );
    }
    
    if (fjInputs.size() == 0) {
        cout << "Error: event with no final state particles" << endl;
        return false;
    }
    return true;
}

void PTCut::JetLoop()
{
    int counter = 0;
    TLorentzVector zero(0,0,0,0);
    for (unsigned int i = 0; i < sortedJets.size(); i++) {
        if (fabs(sortedJets[i].pseudorapidity()) > etaMax) continue;
        if (counter++ == 2) continue;
        
        jetParts = sorted_by_pt(sortedJets[i].constituents());
        if ( jetParts.size() < 2 ) continue;
        
        // Initialize four-momentum sums
        sumPPt = zero; sumPHCAL = zero; sumP3GeV = zero; sumPHi = zero;
        sumCaloPPt = zero; sumCaloPHCAL = zero; sumCaloP3GeV = zero;
        
        sumHistPt = zero; histCh = zero; histNh = zero; histE = zero;
        histMu = zero; histPh = zero;

        // Loop over particles within a jet
        ParticleLoop();
        
        // Fill the histograms:
        HistFill(i);
    }
}


/* TODO: possible to save four-momentum for CHCAL, NHCAL and others into a root file */
void PTCut::ParticleLoop()
{  
    for ( unsigned int j = 0; j != jetParts.size(); ++j ) {
        int absId = jetParts[j].user_index(); 

        int isNHCAL = std::find(neutrals,neutrals+neutrSize,absId)!=(neutrals+neutrSize);
        int isCHCAL = std::find(chargeds,chargeds+chargSize,absId)!=(chargeds+chargSize); 
        
        TLorentzVector tmpP( jetParts[j].px(), jetParts[j].py(), jetParts[j].pz(), jetParts[j].e() );
        double RHCALf = 1;

        // In case of hadrons, set initial values for temporary containers
        if ( isNHCAL || isCHCAL ) {
            RHCALf = fp1->Eval(max(7.0,tmpP.Perp()));  
            RHCAL->Fill( tmpP.Perp(), RHCALf ); // Record fractions
        }

        if ( isCHCAL ) {
            /* Charged hadrons */
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
        } else if ( isNHCAL ) {
            /* Neutral hadrons */
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

            double histResp = HistResp(tmpP, respHCAL, respECAL, respEHHCAL, respEHECAL);
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
        } else {
            /* Others */
            bool known = true;
            if ( absId == 22 ){
                histPh += tmpP;
            } else if ( absId == 11 ){
                histE += tmpP;
            } else if ( absId == 13 ){
                histMu += tmpP;
            } else {
                known = false;
            }
            if (known) sumHistPt += tmpP;
                
            sumP3GeV += tmpP; 
            sumCaloP3GeV += tmpP;
            sumPHCAL += tmpP; 
            sumCaloPHCAL += tmpP;
            sumPPt += tmpP;
            sumCaloPPt += tmpP; 
            sumPHi += tmpP;
        }
    }
}

void PTCut::HistFill(int i)
{
    histProfile->Fill( sortedJets[i].pt(), sumHistPt.Perp()/sortedJets[i].pt() );
    histPNh->Fill( sortedJets[i].pt(), histNh.Perp()/sumHistPt.Perp() );
    histPCh->Fill( sortedJets[i].pt(), histCh.Perp()/sumHistPt.Perp() );
    histPE->Fill( sortedJets[i].pt(), histE.Perp()/sumHistPt.Perp() );
    histPMu->Fill( sortedJets[i].pt(), histMu.Perp()/sumHistPt.Perp() );
    histPPh->Fill( sortedJets[i].pt(), histPh.Perp()/sumHistPt.Perp() );
    
    ptProfile->Fill( sortedJets[i].pt(), sumPPt.Perp()/sortedJets[i].pt() );
    hcalProfile->Fill( sortedJets[i].pt(), sumPHCAL.Perp()/sortedJets[i].pt() );
    gev3Profile->Fill( sortedJets[i].pt(), sumP3GeV.Perp()/sortedJets[i].pt() );
    ptCaloProfile->Fill( sortedJets[i].pt(), sumCaloPPt.Perp()/sortedJets[i].pt() );
    hcalCaloProfile->Fill( sortedJets[i].pt(), sumCaloPHCAL.Perp()/sortedJets[i].pt() );
    gev3CaloProfile->Fill( sortedJets[i].pt(), sumCaloP3GeV.Perp()/sortedJets[i].pt() );
    hiEProfile->Fill( sortedJets[i].pt(), sumPHi.Perp()/sortedJets[i].pt() );
}


void PTCut::WriteResults()
{
    outFile->cd();
    
    ptProfile->Write();
    hcalProfile->Write();
    gev3Profile->Write();
    ptCaloProfile->Write();
    hcalCaloProfile->Write();
    gev3CaloProfile->Write();
    hiEProfile->Write();
    jetMultipl->Write();

    histProfile->Write();
    histPNh->Write();
    histPCh->Write();
    histPE->Write();
    histPMu->Write();
    histPPh->Write();
    
    RHCAL->Write();
    
    outFile->Close();
}


double PTCut::HistResp(TLorentzVector tmpP, double respHCAL, double respECAL,
                         double respEHHCAL, double respEHECAL)
{
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