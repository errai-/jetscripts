// Purpose: Fill inclusive jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: Aug 9, 2011

#define fillHistos_cxx
#include "fillHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TDirectory.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

void fillHistos::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fillHistos.C
//      Root > fillHistos t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t ntot = fChain->GetEntries();
   Long64_t nskip = 0;
   //nentries = 1000;//very short test runs
   //nentries = 100000;//short test runs
   //nentries = 1000000;//medium test runs
   //nskip = 8957000;//8900000;//10000000;//short test runs
   //assert(nentries+nskip <= fChain->GetEntriesFast());
   //assert(nentries+nskip <= ntot);//fChain->GetEntriesFast());

   map<string, int> cnt; // efficiency counters

   //TStopwatch t;
   //t.Start();

   ferr = new ofstream("error.log",ios::out);
     
   // Report memory usage to avoid malloc problems when writing file
   *ferr << endl << "Starting Loop() initialization:" << endl << flush;
   cout << endl << "Starting Loop() initialization:" << endl << flush;
   MemInfo_t info;
   gSystem->GetMemInfo(&info);
   *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
   cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

   _entries = (_mc ? fChain->GetEntries() : ntot);
   _xsecMinBias = 7.126E+10;
   _nbadevts_dup = _nbadevts_run = _nbadevts_ls = _nbadevts_lum = 0;
   _nbadevts_veto = _nbadevts_stream = 0;
   _bscounter_bad = _bscounter_good = _halocounter_bad = _halocounter_good = 0;
   _ecalcounter_good = _ecalcounter_bad = 0;
   _rhocounter_good = _rhocounter_bad = 0;
   _trgcounter = _evtcounter = _totcounter = 0;
   ecalveto = 0;

   // Set cross section weights for pThat bins
   hmcweight = 0;
   if (_pthatbins) {
     double bins[]={5,15,30,50,80,120,170,300,470,600,800,1000,1400,1800,3500};
     const int nbins = sizeof(bins)/sizeof(bins[0])-1;
     double lums[] = {4.49e-5,1.32e-2,1.22e-1,1.04,8.40,5.32e+1,2.57e+2,5.51e+3,
		      5.68e+4,2.73e+5,2.20e+6,6.30e+6,2.02e+8,8.20e+8};
     const int nlums = sizeof(lums)/sizeof(lums[0]);
     assert(nbins==nlums);
     hmcweight = new TH1D("hmcweight",";#hat{p}_{T} (GeV)",nbins,bins);
     for (int i = 0; i != nbins; ++i) {
       hmcweight->SetBinContent(i+1, 1./lums[i]);
     } // for i
   }

   if (_quick) {

     fChain->SetBranchStatus("*",0);

     // Luminosity calculation
     if (_mc) fChain->SetBranchStatus("EvtHdr_.mPthat",1); // pthat
     if (_mc) fChain->SetBranchStatus("EvtHdr_.mWeight",1); // weight
     if (_rd) fChain->SetBranchStatus("EvtHdr_.mRun",1); // run
     if (_rd) fChain->SetBranchStatus("EvtHdr_.mEvent",1); // evt
     if (_rd) fChain->SetBranchStatus("EvtHdr_.mLumi",1); // lbn

     // Event properties
     fChain->SetBranchStatus("EvtHdr_.mNVtx",1); // npv
     fChain->SetBranchStatus("EvtHdr_.mNVtxGood",1); // npvgood
     fChain->SetBranchStatus("EvtHdr_.mPFRho",1); // rho

     // Jet properties (jtpt, jte, jteta, jty, jtphi etc.)
     fChain->SetBranchStatus("PFJets_",1); // njt
     fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fX",1); // jtp4x
     fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fY",1); // jtp4y
     fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fZ",1); // jtp4z
     fChain->SetBranchStatus("PFJets_.P4_.fCoordinates.fT",1); // jtp4t
     fChain->SetBranchStatus("PFJets_.cor_",1); // jtjes
     fChain->SetBranchStatus("PFJets_.area_",1); // jta

     if (_mc) {
       fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fX",1); // jtgenp4x
       fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fY",1); // jtgenp4y
       fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fZ",1); // jtgenp4z
       fChain->SetBranchStatus("PFJets_.genP4_.fCoordinates.fT",1); // jtgenp4t
       fChain->SetBranchStatus("PFJets_.genR_",1); // jtgenr
     }

     // Component fractions
     fChain->SetBranchStatus("PFJets_.chf_",1); // jtchf
     fChain->SetBranchStatus("PFJets_.phf_",1); // jtnef
     fChain->SetBranchStatus("PFJets_.nhf_",1); // jtnhf
     fChain->SetBranchStatus("PFJets_.elf_",1); // jtcef !!
     fChain->SetBranchStatus("PFJets_.muf_",1); // jtmuf !!
     fChain->SetBranchStatus("PFJets_.ncand_",1); // jtn
     fChain->SetBranchStatus("PFJets_.beta_",1); // jtbeta
     fChain->SetBranchStatus("PFJets_.betaStar_",1); // jtbetastar
     fChain->SetBranchStatus("PFJets_.chm_",1); // jtnch
     fChain->SetBranchStatus("PFJets_.phm_",1); // jtnne
     fChain->SetBranchStatus("PFJets_.nhm_",1); // jtnnh
     fChain->SetBranchStatus("PFJets_.elm_",1); // jtnce !!
     fChain->SetBranchStatus("PFJets_.mum_",1); // jtnmu !!
     fChain->SetBranchStatus("PFJets_.tightID_",1); // jtidtight
     fChain->SetBranchStatus("PFJets_.looseID_",1); // jtidloose

     //fChain->SetBranchStatus("rho",1);
     fChain->SetBranchStatus("PFMet_.et_",1); // met
     fChain->SetBranchStatus("PFMet_.phi_",1); // metphi
     fChain->SetBranchStatus("PFMet_.sumEt_",1); // metsumet

     fChain->SetBranchStatus("TriggerDecision_",1);
     fChain->SetBranchStatus("L1Prescale_",1);
     fChain->SetBranchStatus("HLTPrescale_",1);

     // Event cleaning
     //fChain->SetBranchStatus("pvrho",1);
     fChain->SetBranchStatus("EvtHdr_.mPVx",1); // pvx
     fChain->SetBranchStatus("EvtHdr_.mPVy",1); // pvy
     fChain->SetBranchStatus("EvtHdr_.mPVz",1); // pvz
     fChain->SetBranchStatus("EvtHdr_.mPVndof",1); // pvndof
     fChain->SetBranchStatus("EvtHdr_.mBSx",1); // bsx
     fChain->SetBranchStatus("EvtHdr_.mBSy",1); // bsy
     //
     if (_mc) fChain->SetBranchStatus("EvtHdr_.mINTPU",1); // itpu
     if (_mc) fChain->SetBranchStatus("EvtHdr_.mOOTPULate",1); // ootpulate
     if (_mc) fChain->SetBranchStatus("EvtHdr_.mOOTPUEarly",1); // ootpuearly

     if (_mc) {
       fChain->SetBranchStatus("GenJets_",1); // gen_njt
       fChain->SetBranchStatus("GenJets_.fCoordinates.fX",1); // gen_jtp4x
       fChain->SetBranchStatus("GenJets_.fCoordinates.fY",1); // gen_jtp4y
       fChain->SetBranchStatus("GenJets_.fCoordinates.fZ",1); // gen_jtp4z
       fChain->SetBranchStatus("GenJets_.fCoordinates.fT",1); // gen_jtp4t
     }

     if (dofriends) {
       fChain->SetBranchStatus("CaloJets_",1); // njt
       fChain->SetBranchStatus("CaloJets_.P4_.fCoordinates.fX",1); // jtp4x
       fChain->SetBranchStatus("CaloJets_.P4_.fCoordinates.fY",1); // jtp4y
       fChain->SetBranchStatus("CaloJets_.P4_.fCoordinates.fZ",1); // jtp4z
       fChain->SetBranchStatus("CaloJets_.P4_.fCoordinates.fT",1); // jtp4t
       fChain->SetBranchStatus("CaloJets_.cor_",1); // jtjes
       //fChain->SetBranchStatus("CaloJets_.area_",1); // jta
       fChain->SetBranchStatus("CaloJets_.emf_",1); // emf
       fChain->SetBranchStatus("CaloJets_.tightID_",1); // jtidtight
       fChain->SetBranchStatus("CaloJets_.looseID_",1); // jtidloose

       fChain->SetBranchStatus("EvtHdr_.mCaloRho",1); // c_rho

       fChain->SetBranchStatus("CaloMet_.et_",1); // c_met
       fChain->SetBranchStatus("CaloMet_.sumEt_",1); // c_metsumet
     }

     // MC truth information
     //if (_mc) fChain->SetBranchStatus("ak5gen.njt",1);
     //if (_mc) fChain->SetBranchStatus("ak5gen.jtpt",1);
     //if (_mc) fChain->SetBranchStatus("ak5gen.jty",1);
     //if (_mc) fChain->SetBranchStatus("ak5gen.jtgenflv",1);
   } // _quick
   else
     fChain->SetBranchStatus("*",1);

   // Set pointers to branches
   jtp4x = &PFJets__P4__fCoordinates_fX[0];
   jtp4y = &PFJets__P4__fCoordinates_fY[0];
   jtp4z = &PFJets__P4__fCoordinates_fZ[0];
   jtp4t = &PFJets__P4__fCoordinates_fT[0];
   jta = &PFJets__area_[0];
   jtjes = &PFJets__cor_[0];
   jtbeta = &PFJets__beta_[0];
   jtbetastar = &PFJets__betaStar_[0];
   jtidloose = &PFJets__looseID_[0];
   jtidtight = &PFJets__tightID_[0];
   //
   jtgenr = &PFJets__genR_[0];
   jtgenp4x = &PFJets__genP4__fCoordinates_fX[0];
   jtgenp4y = &PFJets__genP4__fCoordinates_fY[0];
   jtgenp4z = &PFJets__genP4__fCoordinates_fZ[0];
   jtgenp4t = &PFJets__genP4__fCoordinates_fT[0];
   //
   c_jtp4x = &CaloJets__P4__fCoordinates_fX[0];
   c_jtp4y = &CaloJets__P4__fCoordinates_fY[0];
   c_jtp4z = &CaloJets__P4__fCoordinates_fZ[0];
   c_jtp4t = &CaloJets__P4__fCoordinates_fT[0];
   //c_jta = &CaloJets__area_[0]; // MISSING!
   //c_jta = &PFJets__area_[0]; // TEMPORARY FIX
   c_jta = TMath::Pi()*0.5*0.5;
   c_jtjes = &CaloJets__cor_[0];
   c_jtemf = &CaloJets__emf_[0];
   c_jtidloose = &CaloJets__looseID_[0];
   c_jtidtight = &CaloJets__tightID_[0];
   //
   jtn = &PFJets__ncand_[0];
   jtnch = &PFJets__chm_[0];
   jtnnh = &PFJets__nhm_[0];
   jtnne = &PFJets__phm_[0];
   jtnce = &PFJets__elm_[0];
   jtnmu = &PFJets__mum_[0];
   jtchf = &PFJets__chf_[0];
   jtnhf = &PFJets__nhf_[0];
   jtnef = &PFJets__phf_[0];
   jtcef = &PFJets__elf_[0];
   jtmuf = &PFJets__muf_[0];
   //
   gen_jtp4x = &GenJets__fCoordinates_fX[0];
   gen_jtp4y = &GenJets__fCoordinates_fY[0];
   gen_jtp4z = &GenJets__fCoordinates_fZ[0];
   gen_jtp4t = &GenJets__fCoordinates_fT[0];

   assert(_algo=="AK5" || _algo=="AK7");
   const char *a = _algo.c_str();
   cout << "\nCONFIGURATION DUMP:" << endl;
   cout << "-------------------" << endl;
   cout << Form("Running over %sPF and %sCALO",a,a) << endl;
   cout << (_useIOV ? "Using" : "Not using")
	<< " time dependent JEC (IOV)" << endl;
   cout << (_doEras ? "Storing all " : "Not storing")
	<< " eras separately" << endl;
   cout << (_doECALveto ? "Vetoing" : "Not vetoing")
	<< " jets in bad ECAL towers" << endl;
   cout << (_doCHS ? "Applying" : "Not applying")
	<< " CHS through betaStar" << endl;
   if (_mc) {
     cout << (_pthatbins ? "Processing pThat binned samples"
	      : "Processing \"flat\" samples") << endl;
   }
   cout << endl;

   // Time dependent JEC
   iov = new jec::IOV(Form("%sPF",a));
   bool isdata = !_mc;
   iov->add("May10",160431,163869,isdata);
   iov->add("PrReV4",165088,167913,isdata);
   iov->add("Aug05",170722,172619,isdata);
   iov->add("PrReV6",172620,173692,isdata);
   iov->add("11BPrV1",175860,180252,isdata);
   iovc = new jec::IOV(Form("%sCalo",a));
   iovc->add("May10",160431,163869,isdata);
   iovc->add("PrReV4",165088,167913,isdata);
   iovc->add("Aug05",170722,172619,isdata);
   iovc->add("PrReV6",172620,173692,isdata);
   iovc->add("11BPrV1",175860,180252,isdata);

   // Full redoing of JEC
   _JEC = 0;
   {
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L1FastJet.txt");
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1FastJet_AK5PF.txt"); // last
     //JetCorrectorParameters *par_l1chs = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1FastJet_AK5PFchs.txt"); // last
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PF_SCALED.txt");
     //JetCorrectorParameters *par_l1chs = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PFchs_SCALED.txt");
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PF_NOMINAL.txt");
     //JetCorrectorParameters *par_l1chs = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PFchs_NOMINAL.txt");
     //JetCorrectorParameters *par_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L2Relative.txt");
     //JetCorrectorParameters *par_l3= new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L3Absolute.txt");
     //JetCorrectorParameters *par_l2l3= new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L2L3Monolithic_AK5PFL1.txt"); // Nov 14 - v0
     //JetCorrectorParameters *par_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_Jec11_V7_AK5PF_L2Relative.txt"); // Nov 28 - last
     //JetCorrectorParameters *par_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_Jec11_V7_AK5PF_L3Absolute.txt"); // Nov 28 - last
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L2L3Residual.txt");
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V3_AK5PF_L2L3Residual.txt"); // JER bias added
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/MPF_METJERcorr_HF_PTDEP_pt_L2L3Residual_AK5PF.txt"); // JER bias + pT dependence added
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/base_RR_PTDEP_pt_L2L3Residual_AK5PF.txt"); // JER bias + pT dependence from pT balance, Nov15
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/L2Res_pT_pTdep.txt"); // JER bias + pT dependence from pT balance, Nov14
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/RR_L2L3Residual_AK5PF.txt"); // JER bias + pT dependence from pT balance, Nov 28
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/11DB_MPF_L2L3Residual_AK5PF.txt"); // Dec 10
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/11DB_MPF_imp_kFSRAbs_AbsEtaFIX_L2L3Residual_AK5PF.txt"); // Dec 12
     // New JEC released on hn-cms-jes on Jan 2nd (GR_R_42_V23)
     JetCorrectorParameters *par_l1 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%sPF.txt",a));
     JetCorrectorParameters *par_l1chs = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%sPFchs.txt",a));
     JetCorrectorParameters *par_l2 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L2Relative_%sPF.txt",a));
     JetCorrectorParameters *par_l3 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L3Absolute_%sPF.txt",a));
     JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L2L3Residual_%sPF.txt",a));
     vector<JetCorrectorParameters> vpar;
     if (_rd && _doCHS) vpar.push_back(*par_l1chs); // doCHS
     else               vpar.push_back(*par_l1);
     vpar.push_back(*par_l2);
     vpar.push_back(*par_l3);
     //vpar.push_back(*par_l2l3);
     if (_rd) vpar.push_back(*par_l2l3res);
     _JEC = new FactorizedJetCorrector(vpar);
     
   } // JEC redone
   assert(_JEC);

   _JEC_ak5calo = 0;
   {
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5Calo_SCALED.txt"); // last
     //JetCorrectorParameters *par_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5Calo_NOMINAL.txt");
     //JetCorrectorParameters *par_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5Calo_L2Relative.txt");
     //JetCorrectorParameters *par_l3= new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5Calo_L3Absolute.txt");
     //JetCorrectorParameters *par_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_Jec11_V7_AK5Calo_L2Relative.txt"); // Nov 28 - last
     //JetCorrectorParameters *par_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_Jec11_V7_AK5Calo_L3Absolute.txt"); // Nov 28 - last
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5Calo_L2L3Residual.txt");
     //
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/MPF_METJERcorr_HF_PTDEP_pt_L2L3Residual_AK5PF.txt"); // JER bias + pT dependence added, Nov14
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/base_RR_PTDEP_pt_L2L3Residual_AK5PF.txt"); // JER bias + pT dependence from pT balance, Nov15
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/RR_L2L3Residual_AK5Calo.txt"); // JER bias + pT dependence from pT balance, Nov 28
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/11DB_MPF_L2L3Residual_AK5Calo.txt"); // Dec 10
     //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/11DB_MPF_imp_kFSRAbs_AbsEtaFIX_L2L3Residual_AK5Calo.txt"); // Dec 12 - last
     // New JEC released on Jan 2dn (GR_R_42_V23)
     JetCorrectorParameters *par_l1 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1Offset_%sCalo.txt",a));
     JetCorrectorParameters *par_l2 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L2Relative_%sCalo.txt",a));
     JetCorrectorParameters *par_l3 = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L3Absolute_%sCalo.txt",a));
     JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L2L3Residual_%sCalo.txt",a));
     vector<JetCorrectorParameters> vpar;
     vpar.push_back(*par_l1);
     vpar.push_back(*par_l2);
     vpar.push_back(*par_l3);
     if (_rd) vpar.push_back(*par_l2l3res);
     _JEC_ak5calo = new FactorizedJetCorrector(vpar);
     
   } // Calo JEC redone
   assert(_JEC_ak5calo);

   // JEC for L1 pile-up studies
   {
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L1Offset.txt");
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PF_SCALED.txt"); // last
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PF_NOMINAL.txt");
     JetCorrectorParameters *L1OffPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1Offset_%sPF.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1OffPar);
     _L1Off = new FactorizedJetCorrector(vParam);
   }
   {
     //JetCorrectorParameters *L1FastPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L1FastJet.txt");
     //JetCorrectorParameters *L1FastPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1FastJet_AK5PF.txt"); // last
     JetCorrectorParameters *L1FastPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%sPF.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1FastPar);
     _L1Fast = new FactorizedJetCorrector(vParam);
   }
   /*
   {
     JetCorrectorParameters *L1FastPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_L1FastJet.txt");
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1FastPar);
     _L1FastOld = new FactorizedJetCorrector(vParam);
   }
   */
   {
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PFchs_SCALED.txt"); // last
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5PFchs_NOMINAL.txt");
     JetCorrectorParameters *L1OffPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1Offset_%sPFchs.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1OffPar);
     _L1OffCHS = new FactorizedJetCorrector(vParam);
   }
   {
     //JetCorrectorParameters *L1FastPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1FastJet_AK5PFchs.txt");
     JetCorrectorParameters *L1FastPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%sPFchs.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1FastPar);
     _L1FastCHS = new FactorizedJetCorrector(vParam);
   }

   {
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5Calo_SCALED.txt"); // last
     //JetCorrectorParameters *L1OffPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1Offset_AK5Calo_NOMINAL.txt");
     JetCorrectorParameters *L1OffPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1Offset_%sCalo.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1OffPar);
     _L1Off_ak5calo = new FactorizedJetCorrector(vParam);
   }
   {
     //JetCorrectorParameters *L1FastPar = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer11_L1FastJet_AK5Calo.txt");
     JetCorrectorParameters *L1FastPar = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%sCalo.txt",a));
     vector<JetCorrectorParameters> vParam;
     vParam.push_back(*L1FastPar);
     _L1Fast_ak5calo = new FactorizedJetCorrector(vParam);
   }

   _jecUnc = 0;
   //_jecUnc2 = 0;
   if (_rd) {
     //jecUnc = new JetCorrectionUncertainty("CondFormats/JetMETObjects/data/JEC11_V2_AK5PF_Uncertainty.txt");
     //_jecUnc = new JetCorrectionUncertainty("CondFormats/JetMETObjects/data/JECUncert2011_AK5PF.txt"); // last
   //_jecUnc2 = new L3Corr(L3Corr::PFAK5_DATA,jec::kAll,L3Corr::kMedium,true);
     _jecUnc = new JetCorrectionUncertainty(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_Uncertainty_%sPF.txt",a));
   }


   // Load latest JSON selection
   if (_rd) {
     //loadJSON("lumicalc/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON.txt");
     //loadJSON("lumicalc/Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON.txt");
     //loadJSON("lumicalc/Cert_160404-166861_7TeV_PromptReco_OR_May10ReReco_JSON.txt");
     //loadJSON("lumicalc/Cert_160404-172255_7TeV_PromptReco_OR_May10ReReco_JSON.txt");
     //loadJSON("lumicalc/Cert_160404-173692_7TeV_May10ReReco_OR_Aug5ReReco_OR_PromptReco_JSON.txt");
     //loadJSON("lumicalc/Oct11th/lumiSummary_Jet_Oct11th.json");
     //loadJSON("lumicalc/Nov13th/lumiSummary.json");
     loadJSON("lumicalc/lumiSummary_Jan16th.json");
     // May10ReReco has few more runs marked as good
   }

   // Load PU profiles for MC reweighing
   //loadPUProfiles("pileup/pudist.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Oct11th.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Oct11th_pertrig.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Oct11th_mb735.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Nov13th_mb735_partial.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Nov13th_mb717.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Nov13th_mb735.root","pileup/pileup_PY.root");
   //loadPUProfiles("pileup/pileup_Jan16th_mb680.root","pileup/pileup_PY.root");
   loadPUProfiles("pileup/pileup_Jan16th_mb735.root","pileup/pileup_PY.root");
   
   // load ECAL veto file for cleaning data
   //if (_rd) {
   loadECALveto("lumicalc/ecalveto.root");
   //}
   
   // Add these runs to the manual veto list
   if (_rd) {

     // Veto list for 38X 36/pb
     //_runveto.insert(142418);
     // _runveto.insert(162765); // rate 30% in all triggers for 1.1/fb
   }

   // load luminosity tables (prescales now stored in event)
   if (_rd) {
     //loadLumi("lumicalc/lumicalc2_by_LS.csv");
     //loadLumi("lumicalc/lumicalc2v1_by_LS.csv"); // Jan16th JSON
     //loadLumi("lumicalc/lumicalc2v2_by_LS.csv"); // Jan16th JSON
     //loadLumi("lumicalc/lumicalc2v3_by_LS.csv"); // Jan16th JSON
     loadLumi("lumicalc/pixellumi_by_LS.csv"); // Jan16th JSON; Feb 22
   }
   if (_mc) cout << Form("Running on MC produced with %1.3g nb-1 (%ld evts)",
			 1000. * _entries / _xsecMinBias,
			 (long int)_entries) << endl;
   if (_rd) cout << Form("Running on %ld events of data",
			 (long int)_entries) << endl;

   // Initialize histograms for different epochs and DQM selections
   initBasics("Standard");
   if (_rd && _doEras) {
     initBasics("May10");
     initBasics("PromptV4");
     initBasics("Aug5");
     initBasics("PromptV6"); // bug: V6 on jan4
     initBasics("2011B");
   }

   if (_rd) {
     initRunHistos("Runs",0.,3.);
     initRunHistos("RunsBarrel",0.,1.);
     initRunHistos("RunsTransition",1.,2.);
     initRunHistos("RunsEndcap",2.,3.);
   }

   // Report memory usage to avoid malloc problems when writing file
   *ferr << "Beginning Loop() proper:" << endl << flush;
   cout  << "Beginning Loop() proper:" << endl << flush;
   gSystem->GetMemInfo(&info);
   *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
   cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

   TStopwatch t;
   t.Start();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=nskip; jentry<nentries+nskip;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%50000==0) cout << "." << flush;

      if (jentry==10000 || jentry==100000 || jentry==1000000 || jentry==5000000){
	cout << endl
	     << Form("Processed %ld events (%1.1f%%) in %1.0f sec. ETA:",
		     (long int)jentry, 100.*jentry/ntot,
		     t.RealTime()) << endl;
	TDatime now; now.Set(now.Convert()+t.RealTime()*ntot/jentry);
	now.Print();
	t.Continue();
      }

      // Set auxiliary event variables (jets, triggers later)
      assert(!_pthatbins || hmcweight);
      pthat = EvtHdr__mPthat;
      weight = (_pthatbins ?
		hmcweight->GetBinContent(hmcweight->FindBin(pthat)) :
		EvtHdr__mWeight);
      run = EvtHdr__mRun;
      evt = EvtHdr__mEvent;
      lbn = EvtHdr__mLumi;
      itpu = EvtHdr__mINTPU;
      ootpulate = EvtHdr__mOOTPULate;
      ootpuearly = EvtHdr__mOOTPUEarly;

      npv = EvtHdr__mNVtx;
      npvgood = EvtHdr__mNVtxGood;
      pvx = EvtHdr__mPVx;
      pvy = EvtHdr__mPVy;
      pvz = EvtHdr__mPVz;
      pvndof = EvtHdr__mPVndof;
      bsx = EvtHdr__mBSx;
      bsy = EvtHdr__mBSy;

      rho = EvtHdr__mPFRho;
      met = PFMet__et_;
      metphi = PFMet__phi_;
      metsumet = PFMet__sumEt_;

      c_rho = EvtHdr__mCaloRho;
      c_met = CaloMet__et_;
      c_metsumet = CaloMet__sumEt_;

      njt = PFJets__;       //assert(njt < kMaxPFJets_);
      c_njt = CaloJets__;   //assert(c_njt < kMaxCaloJets_);
      if (!(njt <= kMaxPFJets_ && c_njt <= kMaxCaloJets_)) {
	*ferr << "Array overflow: njt = "<<njt<<", c_njt = "<<c_njt<<endl;
	cout << "Array overflow: njt = "<<njt<<", c_njt = "<<c_njt<<endl;
	njt = min(njt, PFJets__);
	c_njt = min(c_njt, CaloJets__);
      }
      gen_njt = GenJets__;

      if (_debug) {
	cout << endl << flush;
	Show(jentry);
	cout << endl << flush;
      }

      // Check if duplicate
      if (_rd && _checkduplicates) {
	set<int>& events = _duplicates[run][lbn];
	if (events.find(evt)!=events.end()) {
	  ++_nbadevts_dup;
	  continue;
	}
	events.insert(evt);
      }      
      
      ++cnt["01all"];

      // Check if good run/LS, including JSON selection
      if (_rd) {

	// Does the run/LS pass the latest JSON selection?
	if (_json[run][lbn]==0) {
	  _badjson.insert(make_pair<int, int>(run, lbn));
	  ++_nbadevts_json;
	  if (_dojson) continue;
	}

	// Do we have the run listed in the .csv file?
	map<int, map<int, float> >::const_iterator irun = _lums.find(run);
	if (irun==_lums.end()) {
	  _badruns.insert(run);
	  ++_nbadevts_run;
	  continue;
	}
	// Do we have the LS listed in the .csv file?
	map<int, float>::const_iterator ils = irun->second.find(lbn);
	if (ils==irun->second.end()) {
	  _badlums.insert(make_pair<int, int>(run,lbn));
	  ++_nbadevts_ls;
	  continue;
	}
	// Does the .csv file list a non-zero luminosity?
	if (ils->second==0) {
	  _nolums.insert(make_pair<int, int>(run, lbn));
	  ++_nbadevts_lum;
	  //continue; // Could be Poisson fluctuation to zero
	}

	// Do we exercise run veto based on cross section stability?
	if (_runveto.find(run)!=_runveto.end()) {
	  ++_nbadevts_veto;
	  continue;
	}

	// Keep track of LBNs
	_jt15lums.insert(make_pair<int, int>(run, lbn));
      }

      // Reset event ID
      _pass = true;

      if (_pass) ++cnt["02ls"];

      // Reject events with no vertex
      pvrho = tools::oplus(pvx, pvy);
      _pass = (npvgood>0 && pvrho<2.);

      if (_pass) ++cnt["03vtx"];

      // Event cuts against beam backgrounds
      if (_pass && (tools::oplus(pvx-bsx, pvy-bsy)>0.15 ||
		    pvndof<=4 || fabs(pvz) >= 24.)) {
	++_bscounter_bad;
	_pass = false;
      }
      if (_pass) ++_bscounter_good;
      if (_pass) ++cnt["04bsc"];

      // Event cuts against beam backgrounds
      /* // Commented out on Aug 4 for testing
      if (_pass && ecalveto &&
	  ((njt>=1 &&
	    ecalveto->GetBinContent(ecalveto->FindBin(jteta[0],jtphi[0]))!=0) ||
	   (njt>=2 &&
	    ecalveto->GetBinContent(ecalveto->FindBin(jteta[1],jtphi[1]))!=0))){
	++_ecalcounter_bad;
	_pass = false;
      } // ecal veto
      */
      if (_pass) ++_ecalcounter_good;
      if (_pass) ++cnt["05ecal"];

      // Check rho
      if (_pass && rho>40.) {
	++_rhocounter_bad;
	_pass = false;
	cout << Form("\nrun:ev:ls %d:%d:%d : rho=%1.1f njt=%d npv=%d"
		     " jtpt0=%1.1f sumet=%1.1f met=%1.1f\n",
		     run, lbn, evt, rho, njt, npv,
		     (njt>0 ? jtpt[0] :0.), metsumet, met) << flush;
      }
      if (_pass) ++_rhocounter_good;
      if (_pass) ++cnt["06rho"];


      // Reset prescales (dynamic can change within run)
      for (map<std::string, std::map<int, int> >::iterator it
	     = _prescales.begin(); it != _prescales.end(); ++it) {
	it->second[run] = 0;
      }

      // Fill trigger information
      _trigs.clear();
      if (_mc) _trigs.insert("mc");

      //assert(_mc || TriggerDecision_.size()==35); // 1/fb
      //assert(_mc || TriggerDecision_.size()==53); // 2/fb
      //assert(_mc || TriggerDecision_.size()==54); // 3.2/fb before fix
      //assert(_mc || TriggerDecision_.size()==55); // 3.2/fb after fix
      assert(_mc || TriggerDecision_.size()==64); // 4.7/fb
      for (unsigned int itrg = 0; itrg != TriggerDecision_.size(); ++itrg) {

	bool pass = (TriggerDecision_[itrg]==1); // -1, 0, 1
	string strg = "";
	/* // 1/fb Aug22nd
	if (            itrg<= 3) strg = "jt30";
	if (itrg>= 4 && itrg<= 7) strg = "jt60";
	if (itrg>= 8 && itrg<=11) strg = "jt80";
	if (itrg>=12 && itrg<=15) strg = "jt110";
	if (itrg>=16 && itrg<=19) strg = "jt150";
	if (itrg>=20 && itrg<=23) strg = "jt190";
	if (itrg>=24 && itrg<=27) strg = "jt240";
	if (itrg>=28 && itrg<=30) strg = "jt300"; // v1,2,3 only
	if (itrg>=31 && itrg<=35) strg = "jt370";
	if (itrg>=35) assert(false);
	*/
	/* // 2/fb Sep2nd
	if (            itrg<= 5) strg = "jt30";
	if (itrg>= 6 && itrg<=11) strg = "jt60";
	if (itrg>=12 && itrg<=17) strg = "jt80";
	if (itrg>=18 && itrg<=23) strg = "jt110";
	if (itrg>=24 && itrg<=29) strg = "jt150";
	if (itrg>=30 && itrg<=35) strg = "jt190";
	if (itrg>=36 && itrg<=41) strg = "jt240";
	if (itrg>=42 && itrg<=46) strg = "jt300"; // no v4
	if (itrg>=47 && itrg<=52) strg = "jt370";
	if (itrg>=53) assert(false);
	*/
	/*
	// 3.2/fb Oct11th
	if (            itrg<= 5) strg = "jt30";
	if (itrg>= 6 && itrg<=11) strg = "jt60";
	if (itrg>=12 && itrg<=17) strg = "jt80";
	if (itrg>=18 && itrg<=23) strg = "jt110";
	if (itrg>=24 && itrg<=29) strg = "jt150";
	if (itrg>=30 && itrg<=35) strg = "jt190";
	if (itrg>=36 && itrg<=41) strg = "jt240";
	// Before fix:
	//if (itrg>=42 && itrg<=46) strg = "jt300"; // no v4, add v6 (later)
	//if (itrg==47) {}; // messed up in Oct11th samples?
	//if (itrg>=48 && itrg<=53) strg = "jt370"; // add v7 (lost v1)
	//if (itrg>=54) assert(false);
	// After fix:
	if (itrg>=42 && itrg<=47) strg = "jt300"; // added v6
	if (itrg>=48 && itrg<=54) strg = "jt370"; // added v7
	if (itrg>=55) assert(false);
	*/
	// 4.7/fb Nov13th (jt80 and jt150 removed, 370 to v10, others to v9)
	if (            itrg<= 8) strg = "jt30";
	if (itrg>= 9 && itrg<=17) strg = "jt60";
	if (itrg>=18 && itrg<=26) strg = "jt110";
	if (itrg>=27 && itrg<=35) strg = "jt190";
	if (itrg>=36 && itrg<=44) strg = "jt240";
	if (itrg>=45 && itrg<=53) strg = "jt300";
	if (itrg>=54 && itrg<=63) strg = "jt370"; // up to v10
	if (itrg>=64) assert(false);
	if (pass && strg!="") _trigs.insert(strg);

	// Set prescale from event for now
	if (L1Prescale_[itrg]>0 && HLTPrescale_[itrg]>0)
	  _prescales[strg][run] = max(L1Prescale_[itrg] * HLTPrescale_[itrg],
				      _prescales[strg][run]);

	// check prescale
	if (pass) {
	  double prescale = _prescales[strg][run];
	  if (L1Prescale_[itrg]*HLTPrescale_[itrg]!=prescale) {
	    cout << "Trigger " << strg << ", "
		 << "Prescale(txt file) = " << prescale << endl;
	    cout << "L1 = " << L1Prescale_[itrg] << ", "
		 << "HLT = " << HLTPrescale_[itrg] << endl;
	    assert(false);
	  }
	} // if pass
      } // for itrg

      // Simulate trigger for MC
      //if (_mc && TriggerDecision_.size()==0) {

      //if (c_jtpt[0]>30.) _trigs.insert("jt30");
	// etc.
      //}

      ++_totcounter;
      if (_pass) ++_evtcounter;
      if (_trigs.size()!=0 && _pass) ++_trgcounter;
      if (_trigs.size()!=0 && _pass) ++cnt["07trg"];

      // Retrieve event weight
      _w0 = (_mc ? weight : 1);
      _w = _w0;

      // Calculate trigger weight
      for (unsigned int itrg = 0; itrg != _triggers.size(); ++itrg) {
      
	const char *t = _triggers[itrg].c_str();
	_wt[t] = 1.;

	// Reweigh in-time pile-up
	if (_mc) {
	  int k = pudist[t]->FindBin(itpu);
	  double w1 = pudist[t]->GetBinContent(k);
	  double w2 = pumc->GetBinContent(k);
	  Double_t wpu = (w1==0 || w2==0 ? 1. : w1 / w2);
	  _wt[t] *= wpu;
	}
	// Reweigh late out-of-time pile-up
	if (_mc) {
	  int k = pudist[t]->FindBin(ootpulate);
	  double w1 = pudist[t]->GetBinContent(k);
	  double w2 = pumc->GetBinContent(k);
	  Double_t wpu = (w1==0 || w2==0 ? 1. : w1 / w2);
	  _wt[t] *= wpu;
	}
	// Reweigh early out-of-time pile-up
	if (_mc) {
	  int k = pudist[t]->FindBin(ootpuearly);
	  double w1 = pudist[t]->GetBinContent(k);
	  double w2 = pumc->GetBinContent(k);
	  Double_t wpu = (w1==0 || w2==0 ? 1. : w1 / w2);
	  _wt[t] *= wpu;
	}
      } // for itrg
      _wt["mc"] = _wt["jt370"];

      // To-do: implement reweighing for k-factor (NLO*NP/LOMC)

      // load correct IOV for JEC
      if (_rd && _useIOV) {
	_JEC = iov->get(run); assert(_JEC);
	_JEC_ak5calo = iovc->get(run); assert(_JEC_ak5calo);
      }

      // Calculate pT, eta, phi, y, E and uncorrected pT
      // oversmear jets and MET in MC
      double mex = met * cos(metphi);
      double mey = met * sin(metphi);
      for (int i = 0; i != njt; ++i) {

	// Kostas stores UNCORRECTED four-vector
	p4.SetPxPyPzE(jtp4x[i],jtp4y[i],jtp4z[i],jtp4t[i]);
	jtptu[i] = p4.Pt();
	jteu[i] = p4.E();

	// Apply pseudo-CHS
	//double k = (_mc ? 1 : 1. - jtchf[i]*(1.-jtbeta[i]));
	//double k = ((_mc || !_doCHS) ? 1 : 1. - jtchf[i]*jtbetastar[i]);
	double k = (!_doCHS ? 1 : 1. - jtchf[i]*jtbetastar[i]);
	if (k==0) k = 1e-4; // remove these with jetid
	p4.SetPxPyPzE(k*jtp4x[i], k*jtp4y[i], k*jtp4z[i], k*jtp4t[i]);
	jtptuchs[i] = p4.Pt();
	jteuchs[i] = p4.E();
	//jtchf[i] = jtchf[i]*jtbeta[i]/k;
	if (!_mc && _doCHS) {
	  jtchf[i] = jtchf[i]*(1-jtbetastar[i])/k;
	  jtnef[i] /= k;
	  jtnhf[i] /= k;
	  jtcef[i] /= k;
	  jtmuf[i] /= k;
	  //jtptuchs[i] = jtptu[i]; // doCHS off
	}

	// Recalculate JEC
	_JEC->setRho(rho);
	_JEC->setNPV(npvgood);
	_JEC->setJetA(jta[i]);
	_JEC->setJetPt(jtptuchs[i]);//p4.Pt());
	_JEC->setJetE(jteuchs[i]);// Nov8
	_JEC->setJetEta(p4.Eta());
	//_JEC->setJetEta(max(-4.7, min(4.7, p4.Eta()))); // fix HF JEC
	//jtjesnew[i] = _JEC->getCorrection();
	//jtjesnew[i] = max(_JEC->getCorrection(),(float)1e-4); // Nov 21b
	jtjesnew[i] = _JEC->getCorrection(); // Dec 12
	//jtjesnew[i] *= (_rd ? 0.99 : 1.00);//1.01; // L3 0.985->0.995
	// above now set to 1.01 by default (Nov 14)
	//assert(fabs(jtjesnew[i]-jtjes[i])<1e-4);

	// Recalculate JEC (again to get subcorrections)
	_JEC->setRho(rho);
	_JEC->setNPV(npvgood);
	_JEC->setJetA(jta[i]);
	_JEC->setJetPt(jtptuchs[i]);
	_JEC->setJetE(jteuchs[i]);
	_JEC->setJetEta(p4.Eta());
	//
	vector<float> v = _JEC->getSubCorrections();
	//cout << "v.size() = " << v.size() << endl << flush;
	assert((_rd && v.size()==4) || (_mc && v.size()==3));
	double jec_l1 = v[0];
	double jec_l2l3 = v[2]/v[0];
	double jec_res = (_rd ? v[3]/v[2] : 1.); 
	jtjes_l1[i] = jec_l1;
	jtjes_l2l3[i] = jec_l2l3;
	jtjes_res[i] = jec_res;
	assert(jtjesnew[i] = v[v.size()-1]);
	
	// Calculate corrected compositions
	_L1FastCHS->setRho(rho);
	_L1FastCHS->setJetA(jta[i]);
	_L1FastCHS->setJetPt(jtptu[i]);
	_L1FastCHS->setJetE(jteu[i]);
	_L1FastCHS->setJetEta(p4.Eta());
	double l1 = _L1FastCHS->getCorrection();
	double l1fast = (1 - l1) * jtptu[i];
	double pt_ch = jtptu[i] * jtchf[i] * (1 - jtbetastar[i]);
	double pt_ne = max(jtptu[i] * jtnef[i] - l1fast, 0.);
	double pt_nh = jtptu[i] * jtnhf[i];
	double pt_ce = jtptu[i] * jtcef[i];
	double pt_mu = jtptu[i] * jtmuf[i];
	double jtpt2 = max(pt_ch + pt_ne + pt_nh + pt_ce + pt_mu, 1e-3);
	jtchf2[i] = pt_ch / jtpt2;
	jtnef2[i] = pt_ne / jtpt2;
	jtnhf2[i] = pt_nh / jtpt2;
	jtcef2[i] = pt_ce / jtpt2;
	jtmuf2[i] = pt_mu / jtpt2;

	// Correct jets
	p4 *= jtjesnew[i];
	jte[i] = p4.E();
	jtpt[i] = p4.Pt();
	jteta[i] = p4.Eta();
	jtphi[i] = p4.Phi();
	jty[i] = p4.Rapidity();

	// Calculate gen level info
	if (_mc) {
	  gp4.SetPxPyPzE(jtgenp4x[i],jtgenp4y[i],jtgenp4z[i],jtgenp4t[i]);
	  jtgenpt[i] = gp4.Pt();
	  //jtgeneta[i] = gp4.Eta();//Pt();
	  jtgeny[i] = gp4.Rapidity();//Pt();
	}
	
	// Oversmear MC to match data
	// Results from Matthias Schroeder, JA July 21, 2011:
	// https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=148123
	/* // Turned off on Oct 4 for testing
	if (_mc && jtgenr[i] < 0.25 && jtgenpt[i] > 0. &&
	    jtpt[i] > 0.5*jtgenpt[i] && jtpt[i] < 2.0*jtgenpt[i]) {
	  double kover = 1.00;
	  double x = fabs(jteta[i]);
	  // updated to noResJEC on Aug29 (for |eta|>2.3 in particular)
	  if (x < 0.5)              kover = 1.062;//1.052;
	  if (x >= 0.5 && x < 1.1)  kover = 1.057;//1.057;
	  if (x >= 1.1 && x < 1.7)  kover = 1.089;//1.096;
	  if (x >= 1.7 && x < 2.3)  kover = 1.127;//1.134;
	  if (x >= 2.3 && x < 5.0)  kover = 1.158;//1.288;
	  double dpt = (kover - 1) * (jtpt[i] - jtgenpt[i]);
	  double kf = 1 + dpt / jtpt[i];
	  jte[i] *= kf;
	  jtpt[i] *= kf;
	  jtptu[i] *= kf;
	  // To-do: propagate this to MET also
	  mex -= dpt * cos(jtphi[i]) / jtjesnew[i];
	  //mey -= dpt * cos(jtphi[i]) / jtjesnew[i]; // BUG!!!
	  mey -= dpt * sin(jtphi[i]) / jtjesnew[i]; // BUG!!!
	}
	*/

	met = tools::oplus(mex, mey);
	metphi = atan2(mey, mex);
      } // for i

      if (_mc) {
	for (int i = 0; i != gen_njt; ++i) {

	  genp4.SetPxPyPzE(gen_jtp4x[i],gen_jtp4y[i],gen_jtp4z[i],gen_jtp4t[i]);
	  gen_jtpt[i] = genp4.Pt();
	  gen_jteta[i] = genp4.Eta(); // for matching
	  gen_jtphi[i] = genp4.Phi(); // for matching
	  gen_jty[i] = genp4.Rapidity();
	  // for matching
	} // for i
      } // _mc

      // Propagate jec to MET
      double ucx = -mex;
      double ucy = -mey;
      for (unsigned int i = 0; i != njt; ++i) {
	
	// Only use jets with corr. pT>25 GeV to equalize data and MC thresholds
	if (jtpt[i] > _recopt && fabs(jteta[i])<4.7) {

          // Subtract uncorrected jet pT from met, put back corrected  
          // Remember that MET is negative vector sum
	  // and that uncorrected pT has offset than corrected doesn't
	  _L1Fast->setRho(rho);
	  _L1Fast->setJetA(jta[i]);
	  _L1Fast->setJetPt(jtptu[i]);
	  _L1Fast->setJetE(jteu[i]);
	  _L1Fast->setJetEta(jteta[i]);
	  //_L1Fast->setJetEta(max(-4.7, min(4.7, (double)jteta[i])));
	  double l1corr = _L1Fast->getCorrection();
          double dpt = jtpt[i] - l1corr*jtptu[i];
          mex -= dpt * cos(jtphi[i]);
          mey -= dpt * sin(jtphi[i]);
          // Keep track of remaining pT in unclustered energy, i.e. 
          // subtract jets from -MET to have the non-jet component
	  // treat UE and PU underneath jets as unclustered in order
	  // to keep the homogeneous
	  //double pu = (1 - l1corr) * jtptu[i];
	  //double ue = 1.068/rho * pu;
	  double ue = 1.068 * jta[i];
          ucx -= (l1corr * jtptu[i] - ue) * cos(jtphi[i]);
          ucy -= (l1corr * jtptu[i] - ue) * sin(jtphi[i]);
        }
      } // for i                                                                
      // Type I MET                                                             
      met1 = tools::oplus(mex, mey);
      metphi1 = atan2(mey, mex);
      // Correct unclustered energy; jec for 10 GeV jets varies between         
      // 1.1-1.22 at |y|<2.5, 2.5-3.0 even goes up to 1.35                      
      // => assume average correction of about 1.15 needed                      
      // => did not seem even nearly enough; try 1.5                            
      // => reduce down to 1.25 (high pT threshold on jets)
      mex -= 0.25*ucx;//0.5*ucx;                                             
      //mey -= 0.5*ucx; // BUG!!!
      mey -= 0.25*ucy;//0.5*ucy;
      // Type II MET                                                            
      met2 = tools::oplus(mex, mey);
      metphi2 = atan2(mey, mex);

      // Repeat for CaloJets
      for (int i = 0; i != njt; ++i) {

	// Kostas stores UNCORRECTED four-vector
	cp4.SetPxPyPzE(c_jtp4x[i],c_jtp4y[i],c_jtp4z[i],c_jtp4t[i]);
	c_jtptu[i] = cp4.Pt();
	c_jteu[i] = cp4.E();

	// Recalculate JEC
	//_JEC_ak5calo->setRho(rho); // Bug, Feb 18
	_JEC_ak5calo->setRho(c_rho);
	_JEC_ak5calo->setNPV(npvgood);
	_JEC_ak5calo->setJetA(c_jta);//[i]);
	_JEC_ak5calo->setJetPt(c_jtptu[i]);//cp4.Pt());
	_JEC_ak5calo->setJetE(c_jteu[i]);
	_JEC_ak5calo->setJetEta(cp4.Eta());
	//_JEC_ak5calo->setJetEta(max(-4.7, min(4.7, p4.Eta())));
	//c_jtjesnew[i] = _JEC_ak5calo->getCorrection();
	//c_jtjesnew[i] = max(_JEC_ak5calo->getCorrection(),(float)1e-4);
	c_jtjesnew[i] = _JEC_ak5calo->getCorrection(); // Dec12
	//c_jtjesnew[i] *= (_rd ? 0.99 : 1.00);//1.01; // L3 0.985->0.995
	// above now set to 1.01 by default (Nov 14)

	// Recalculate JEC (again)
	//_JEC_ak5calo->setRho(rho); // Bug, Feb 18
	_JEC_ak5calo->setRho(c_rho);
	_JEC_ak5calo->setNPV(npvgood);
	_JEC_ak5calo->setJetA(c_jta);
	_JEC_ak5calo->setJetPt(c_jtptu[i]);
	_JEC_ak5calo->setJetE(c_jteu[i]);
	_JEC_ak5calo->setJetEta(cp4.Eta());
	vector<float> v = _JEC_ak5calo->getSubCorrections();
	//cout << "vc.size() = " << v.size() << endl << flush;
	assert((_rd && v.size()==4) || (_mc && v.size()==3));
	double jec_l1 = v[0];
	double jec_l2l3 = v[2]/v[0];
	double jec_res = (_rd ? v[3]/v[2] : 1.); 
	c_jtjes_l1[i] = jec_l1;
	c_jtjes_l2l3[i] = jec_l2l3;
	c_jtjes_res[i] = jec_res;
	assert(c_jtjesnew[i] = v[v.size()-1]);

	// Correct jets
	cp4 *= c_jtjesnew[i];
	c_jte[i] = cp4.E();
	c_jtpt[i] = cp4.Pt();
	c_jteta[i] = cp4.Eta();
	c_jtphi[i] = cp4.Phi();
	c_jty[i] = cp4.Rapidity();
      } // for i


      // Fill simulated triggers for MC
      if (_mc && dofriends && TriggerDecision_.size()!=53 && c_njt!=0) {
	if (c_jtpt[0]>30.) _trigs.insert("jt30");
	if (c_jtpt[0]>60.) 
	  _trigs.insert(string("jt60"));
	//_prescales["jt60"][run] = 1;
	//if (c_jtpt[0]>80.) _trigs.insert("jt80");
	if (c_jtpt[0]>110.) _trigs.insert("jt110");
	//if (c_jtpt[0]>150.) _trigs.insert("jt150");
	if (c_jtpt[0]>190.) _trigs.insert("jt190");
	if (c_jtpt[0]>240.) _trigs.insert("jt240");
	if (c_jtpt[0]>300.) _trigs.insert("jt300");
	if (c_jtpt[0]>370.) _trigs.insert("jt370");
      }


      if (njt!=0 && _pass) ++cnt["08njt"];
	
      _jetids.resize(njt);
      for (unsigned int i = 0; i != _jetids.size(); ++i) _jetids[i] = true;
      fillJetID(_jetids);

      if (njt!=0 && _jetids[0] && _pass) ++cnt["09jtid"];

      // Check if overweight PU event
      if (_mc && njt!=0 && _jetids[0] && _pass) {
	//_pass = (jtpt[0] < 1.5*jtgenpt[0] && jtpt[0] < 1.5*pthat);
	_pass = (jtpt[0] < 1.5*jtgenpt[0] || jtpt[0] < 1.5*pthat); // Oct5
	if (_pass) ++cnt["10puw"];
      }

      // Here can categorize events into different triggers, epochs,
      // topologies etc.
      // Eta and pT binning are handled in the fillBasic class 
      fillBasics("Standard");
      if (_rd && _doEras) {
	if (run>=160431 && run<=163869) fillBasics("May10");
	if (run>=165088 && run<=167913) fillBasics("PromptV4");
	if (run>=170722 && run<=172619) fillBasics("Aug5");
	if (run>=172620 && run<=173692) fillBasics("PromptV6");
	//if (run>=175860 && run<=177452) fillBasics("2011B");
	if (run>=175860 && run<=180252) fillBasics("2011B");
      }

      // Run quality checks
      if (_rd) {
	fillRunHistos("Runs");
	fillRunHistos("RunsBarrel");
	fillRunHistos("RunsTransition");
	fillRunHistos("RunsEndcap");
      }

      // Report memory usage to avoid malloc problems when writing file
      if (jentry%1000000==0) {
	*ferr << Form("Doing Loop(), %dM events:",
		      int(jentry/1e6 + 0.5)) << endl << flush;
	gSystem->GetMemInfo(&info);
	*ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d,"
		      " Stot:%d, SUsed:%d, SFree:%d",
		      info.fMemTotal, info.fMemUsed, info.fMemFree,
		      info.fSwapTotal, info.fSwapUsed, info.fSwapFree)
	      << endl << flush;
      } // 1M report

   } // for jentry
   cout << endl;

   // Report memory usage to avoid malloc problems when writing file
   *ferr << "Finished processing " << nentries << " entries:" << endl << flush;
   cout  << "Finished processing " << nentries << " entries:" << endl << flush;
   gSystem->GetMemInfo(&info);
   *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
   cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
		info.fMemTotal, info.fMemUsed, info.fMemFree,
		info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

   writeRunHistos(); // before writeBasics!
   writeBasics();

   // List bad runs
   cout << "Processed " << _totcounter << " events in total" << endl;
   cout << "Processed " << _trgcounter << " events passing "
	<< " basic data quality and trigger cuts" << endl;
   cout << "(out of " << _evtcounter << " passing data quality cuts)" << endl;
   if (_badruns.size()!=0 || _badlums.size()!=0 || _nolums.size()!=0 ||
       _nbadevts_dup!=0 || _nbadevts_json!=0) {
     cout << "Found " << _badruns.size() << " bad runs:";
     for (set<int>::const_iterator it = _badruns.begin();
	  it != _badruns.end(); ++it) {
       cout << " " << *it;
     } // for it
     cout << endl;
     cout << "These contained " << _nbadevts_run << " bad events" << endl;
     cout << "Found " << _nbadevts_json << " bad events according to new JSON"
	  << (_dojson ? " (events cut)" : "(events not cut)") << endl;
     cout << "Found " << _badlums.size() << " bad LS and "
	  << _nolums.size() << " non-normalizable LS in good runs" << endl;
     cout << "These contained " << _nbadevts_ls << " discarded events"
	  << " in bad LS and " << _nbadevts_lum << " in non-normalizable LS"
	  << endl;
     cout << endl;
     cout << "Found " << _nbadevts_dup << " duplicate events, which were"
	  << " properly discarded" << endl;
     cout << "The vetoed runs contained " << _nbadevts_veto
	  << " events" << endl;
   } // has badruns
   cout << "Runs not in JetMETTau stream contained " << _nbadevts_stream
	<< " events" << endl;

   // Report beam spot cut efficiency
   cout << "Beam spot counter discarded " << _bscounter_bad
	<< " events out of " << _bscounter_good
	<< " (" << double(_bscounter_bad)/double(_bscounter_good)*100.
	<< "%)" << endl;
   cout << "Beam spot expectation is less than 0.5%" << endl;

   // Report ECAL hole veto efficiency
   cout << "ECAL hole veto counter discarded " << _ecalcounter_bad
	<< " events out of " << _ecalcounter_good
	<< " (" << double(_ecalcounter_bad)/double(_ecalcounter_good)*100.
	<< "%)" << endl;
   cout << "ECAL hole expectation is less than 2.6% [=2*57/(60*72)]" << endl;

   // Report rho veto efficiency
   cout << "Rho<40 veto counter discarded " << _rhocounter_bad
	<< " events out of " << _rhocounter_good
	<< " (" << double(_rhocounter_bad)/double(_rhocounter_good)*100.
	<< "%)" << endl;
   cout << "Rho veto expectation is less than 1 ppm" << endl;

   // Report beam halo efficiency
   cout << "Beam halo counter flagged (not discarded)" << _halocounter_bad
	<< " events out of " << _halocounter_good
	<< " (" << double(_halocounter_bad)/double(_halocounter_good)*100.
	<< "%) " << endl;
   cout << "This is after the beam spot constraint" << endl;

   cout << endl;
   for (map<string,int>::const_iterator it = cnt.begin(); it != cnt.end(); ++it)
     cout << Form("%s: %d (%1.1f%%)", it->first.c_str(), it->second,
		  100. * it->second / max(1, cnt["01all"])) << endl;
   cout << endl;

   // Report LS actually used for Jet15U in the analysis
   // (not necessarily containing any Jet15U triggers, though)
   cout << "Reporting JetMETTau LS in fillhistos.json" << endl;
   ofstream fout("fillhistos.json", ios::out);
   for (set<pair<int, int> >::const_iterator it = _jt15lums.begin();
	it != _jt15lums.end(); ++it) {
     fout << it->first << " " << it->second << endl;
   }
   if (_dojson) {
     cout << "Reporting LS marked newly bad in fillhistos.json.bad" << endl;
     ofstream fout2("fillhistos.json.bad", ios::out);
     for (set<pair<int, int> >::const_iterator it = _badjson.begin();
	  it != _badjson.end(); ++it) {
       fout2 << it->first << " " << it->second << endl;
     }
   } // _dojson

   t.Stop();
   cout << "Processing used " << t.CpuTime() << "s CPU time ("
	<< t.CpuTime()/3600. << "h)" << endl;
   cout << "Processing used " << t.RealTime() << "s real time ("
	<< t.RealTime()/3600. << "h)" << endl;
   cout << endl << endl;

   delete ferr;
}


// Initialize basic histograms for trigger and eta bins
void fillHistos::initBasics(string name) {

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initBasics("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile :
	      new TFile(Form("output-%s-1.root",_type.c_str()), "RECREATE"));
  assert(f && !f->IsZombie());
  assert(f->mkdir(name.c_str()));
  assert(f->cd(name.c_str()));
  TDirectory *topdir = gDirectory;

  // Rapidity bins + HF + b-tag
  //double y[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 2.0, 2.2};
  // Rapidity bins + HF + barrel
  double y[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 0., 1.3};
  const int ny = sizeof(y)/sizeof(y[0])-1;

  // define triggers
  vector<string> triggers;
  if (_mc) triggers.push_back("mc");
  triggers.push_back("jt30");
  triggers.push_back("jt60");
  //triggers.push_back("jt80");
  triggers.push_back("jt110");
  //triggers.push_back("jt150");
  triggers.push_back("jt190");
  triggers.push_back("jt240");
  triggers.push_back("jt300");
  triggers.push_back("jt370");

  // define efficient pT ranges for triggers for control plots
  map<string, pair<double, double> > pt;
  pt["mc"] = make_pair<double, double>(49., 2000.);
  pt["jt30"] = make_pair<double, double>(56., 84.);
  pt["jt60"] = make_pair<double, double>(84., 153);//114.);
  //pt["jt80"] = make_pair<double, double>(114., 153.);
  pt["jt110"] = make_pair<double, double>(153., 245.);//196.);
  //pt["jt150"] = make_pair<double, double>(196., 245.);
  pt["jt190"] = make_pair<double, double>(245.,300.);
  pt["jt240"] = make_pair<double, double>(300.,395.);
  pt["jt300"] = make_pair<double, double>(395., 468.);
  pt["jt370"] = make_pair<double, double>(468., 2000.);

  map<string, double> pttrg;
  pttrg["mc"] = 20.;
  pttrg["jt30"] = 30.;
  pttrg["jt60"] = 60.;
  //pttrg["jt80"] = 80.;
  pttrg["jt110"] = 110.;
  //pttrg["jt150"] = 150.;
  pttrg["jt190"] = 190.;
  pttrg["jt240"] = 240.;
  pttrg["jt300"] = 300.;
  pttrg["jt370"] = 370.;

  // Loop over rapidity, trigger bins
  for (int i = 0; i != ny; ++i) {

    if (y[i+1] > y[i]) { // create real bins only

      // subdirectory for rapidity bin
      const char *yname = Form("Eta_%1.1f-%1.1f", y[i], y[i+1]);
      assert(topdir);
      assert(topdir->mkdir(yname));
      assert(topdir->cd(yname));
      TDirectory *ydir = gDirectory;
      
      for (unsigned int j = 0; j != triggers.size(); ++j) {
	
	// subdirectory for trigger
	const char *trg = triggers[j].c_str();
	assert(ydir);
	assert(ydir->mkdir(trg));
	assert(ydir->cd(trg));
	TDirectory *dir = gDirectory;

	// Initialize and store
	assert(dir);
	basicHistos *h = new basicHistos(dir, trg, "", y[i], y[i+1], pttrg[trg],
					 pt[trg].first, pt[trg].second,
					 triggers[j]=="mc", dofriends);
	_histos[name].push_back(h);
      } // for j
    } // real bin
  } // for i

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initBasics("<<name<<") finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initBasic

// Loop over basic histogram containers to fill all
void fillHistos::fillBasics(string name) {

  for (unsigned int i = 0; i != _histos[name].size(); ++i)
    fillBasic(_histos[name][i]);
}

// Fill basic histograms after applying pt, y cuts
void fillHistos::fillBasic(basicHistos *h) {

  assert(h);
  h->hpttmp->Reset();

  _w = _w0 * _wt[h->trigname]; assert(_w);

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  // Luminosity information
  if (_rd && h->lums[run][lbn]==0) {
    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    if (lum==0) {
      // cerr << "Run " << run << " LS " << lbn << " empty!" << endl;
      // apparently these can be Poisson fluctuations in the .csv file
    }

    double prescale(0);
    map<int, int>::const_iterator ip = _prescales[h->trigname].find(run);
    if (ip==_prescales[h->trigname].end()) {
      if (fired) {
	*ferr << "No prescale info for trigger " << h->trigname
	      << " in run " << run << "!" << endl << flush; 
	assert(false);
      }
    }
    else
      prescale = ip->second;//_prescales[h->trigname][run];

    if (prescale==0 && fired) {
      *ferr << "Prescale zero for trigger " << h->trigname
	    << " in run " << run << "!" << endl << flush; 
      prescale = 1.;
      assert(false);
    }

    h->lumsum += (prescale ? lum / prescale : 0.);
    h->lumsum2 += (prescale ? lum2 / prescale : 0.);
    //h->lumsum_bts += (prescale && run <= 144114 ? lum / prescale : 0.);
    //h->lumsum_ats += (prescale && run > 144114 ? lum / prescale : 0.);
    h->lums[run][lbn] = 1;
  }

  if (_debug) {
    //if ( ( _mc && string(h->trigname)==string("mc")) ||
    // (!_mc && string(h->trigname)==string("jt30"))) {
    if (h == _histos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (set<string>::iterator it = _trigs.begin();
	   it != _trigs.end(); ++it) {
	cout << *it << ", ";
      }
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired
  //bool fired = (_trigs.find(h->trigname)!=_trigs.end());
  if (!fired) return;
  
  // Check if event looks good
  //bool evtid = (met < 0.5 * metsumet || met < 100.);
  bool evtid = (met < 0.3 * metsumet || met < 50.);

  if (_debug) cout << Form("Subdirectory Eta_%1.1f-%1.1f/%s",
			   h->ymin,h->ymax,h->trigname.c_str()) << endl;
  if (_debug) cout << "Calculate and fill dijet mass" << endl << flush;

  if (h->ismc) h->hpthat->Fill(pthat, _w);
  if (h->ismc) h->hpthatnlo->Fill(pthat);

  if (njt>=2) { // Calculate and fill dijet mass

    // Find leading jets (residual JEC may change ordering)
    map<double, int> ptorder;
    for (int i = 0; i != njt; ++i) {
      double idx = -jtpt[i]; // note minus
      //if (ptorder.find(idx) != ptorder.end()) idx += 1e-5*jteta[i];
      // Nov 21: may have multiple pT=0 jets after L1
      while (ptorder.find(idx) != ptorder.end()) idx += 1e-5*jteta[i];
      assert(ptorder.find(idx)==ptorder.end());
      ptorder[idx] = i;
    }
    int i0 = (ptorder.begin())->second;
    int i1 = (++ptorder.begin())->second;
    //if (i0>1 || i1>1) {
    //cout << "Leading jet pTs in old order:"
    //   << " pt1 = " << _vjec[0]*jtpt[0]
    //   << " pt2 = " << _vjec[1]*jtpt[1]
    //   << " pt3 = " << _vjec[2]*jtpt[2] << endl;
    //cout << "Leading ones now " << i0 << " and " << i1 << endl;
    //} 
    
    //TLorentzVector j1, j2; // use class variables, save instantiation
    j1.SetPtEtaPhiE(jtpt[i0],jteta[i0],jtphi[i0],jte[i0]);
    j2.SetPtEtaPhiE(jtpt[i1],jteta[i1],jtphi[i1],jte[i1]);
    double djmass = (j1+j2).M();
    double ymaxdj = max(fabs(jty[i0]),fabs(jty[i1]));
    bool goodmass = (jtpt[i0]>30. && jtpt[i1]>30.);
    if (evtid && goodmass && _jetids[i0] && _jetids[i1] &&
	ymaxdj >= h->ymin && ymaxdj < h->ymax) {
      assert(h->hdjmass);
      h->hdjmass->Fill(djmass, _w);
      assert(h->hdjmass0);
      h->hdjmass0->Fill(djmass, _w);
      assert(h->pdjmass_ptratio);
      h->pdjmass_ptratio->Fill(djmass, j1.Pt()/j2.Pt(), _w);
      assert(h->pdjmass0_ptratio);
      h->pdjmass0_ptratio->Fill(djmass, j1.Pt()/j2.Pt(), _w);
    }
      
    // Dijet mass for Hgg
    //assert(jtpt[i0]>=jtpt[i1]);
    if (!(jtpt[i0]>=jtpt[i1])) {
      cout << " i0="<<i0<<" pt0="<<jtpt[i0]
	   << " i1="<<i1<<" pt1="<<jtpt[i1] << endl << flush;
    }
    else {
      bool goodhgg = (jtpt[i0]>=40. && jtpt[i1]>=25.);
      if (evtid && goodhgg && _jetids[i0] && _jetids[i1] &&
	  ymaxdj >= h->ymin && ymaxdj < h->ymax) {
	assert(h->hdjmass0_hgg);
	h->hdjmass0_hgg->Fill(djmass, _w);
      }
    }
  } // dijet mass

  if (_debug) cout << "Calculate and fill dijet balance" << endl << flush;

  // Calculate and fill dijet balance histograms
  if (njt>=2 && jtpt[0]>0 && jtpt[1]>0) {
    // Could reorder in pT, but nah (for now)
    int i0 = 0;
    int i1 = 1;
    int i2 = (njt>=3 ? 2 : -1);
    int iref = (fabs(jty[i0]) < fabs(jty[i1]) ? i0 : i1);
    int iprobe = (iref==i0 ? i1 : i0);
    double yref = fabs(jty[iref]);
    double yprobe = fabs(jty[iprobe]);
    double pt1 = jtpt[iref];
    double pt2 = jtpt[iprobe];
    double pt3 = (i2>-1 ? jtpt[i2] : 0.);
    double ptave = 0.5 * (pt2 + pt1); assert(ptave);
    pt3 = (pt3>_recopt ? pt3 / ptave : 0.);
    double asymm = (pt2 - pt1) / (pt2 + pt1);
    double mpf = 1 + met2*cos(delta_phi(metphi2,jtphi[iref]))/ptave;
    //
    double pt3tp = (pt3>_recopt ? pt3 / pt1 : 0.);
    double asymmtp = (pt2 - pt1) / pt1;
    double mpftp = 1 + met2*cos(delta_phi(metphi2,jtphi[iref])) / pt1;

    if (evtid && _jetids[iref] && _jetids[iprobe] && pt1>_recopt && pt2>_recopt
	&& delta_phi(jtphi[iref],jtphi[iprobe]) > 2.8 &&
	yprobe >= h->ymin && yprobe < h->ymax && yref < 1.3) {
      assert(h->hdjasymm);
      h->hdjasymm->Fill(ptave, pt3, asymm, _w);
      assert(h->hdjmpf);
      h->hdjmpf->Fill(ptave, pt3, mpf-1, _w);
      //
      assert(h->hdjasymmtp);
      h->hdjasymmtp->Fill(pt1, pt3tp, asymmtp, _w);
      assert(h->hdjmpftp);
      h->hdjmpftp->Fill(pt1, pt3tp, mpftp-1, _w);
    } // first combo
    // If it so happens that both jets at |y|<1.3, fill also the other combo
    if (evtid && _jetids[iref] && _jetids[iprobe] && pt1>_recopt && pt2>_recopt
	&& delta_phi(jtphi[iref],jtphi[iprobe]) > 2.8 &&
	yref >= h->ymin && yref < h->ymax && yprobe < 1.3) {
      double mpf2 = 1 + met2*cos(delta_phi(metphi2,jtphi[iprobe]))/ptave;
      double asymm2tp = (pt1 - pt2) / pt2;
      double mpf2tp = 1 + met2*cos(delta_phi(metphi2,jtphi[iprobe])) / pt2;
      assert(h->hdjasymm);
      h->hdjasymm->Fill(ptave, pt3, -asymm, _w);
      assert(h->hdjmpf);
      h->hdjmpf->Fill(ptave, pt3, mpf2-1, _w);
      //
      assert(h->hdjasymmtp);
      h->hdjasymmtp->Fill(pt2, pt3tp, asymm2tp, _w);
      assert(h->hdjmpftp);
      h->hdjmpftp->Fill(pt2, pt3tp, mpf2tp-1, _w);
    } // second combo
  }

  // Fill jet pT ratios vs nvtxgood (pile-up)
  bool has2 = (njt>=2 && jtpt[1] > _recopt &&
	       fabs(jty[1])>=h->ymin && fabs(jty[1])<h->ymax);
  bool has3 = (njt>=3 && jtpt[2] > _recopt &&
	       fabs(jty[2])>=h->ymin && fabs(jty[2])<h->ymax
	       && jtpt[1] > 0.70 * jtpt[0]);
  bool has32 = (has3 && fabs(jty[1]) < 1.3);
  if (_pass && evtid && _jetids[0] && jtpt[0]>=h->ptmin && jtpt[0]<h->ptmax &&
      fabs(jty[0]) < 1.3) {
    
    //h->hr21->Fill(npvgood, has2 ? jtpt[1] / jtpt[0] : 0.);
    //h->hr31->Fill(npvgood, has3 ? jtpt[2] / jtpt[0] : 0.);
    //h->hr32->Fill(npvgood, has3 ? jtpt[2] / jtpt[1] : 0.);
    h->hr21->Fill(has2 ? jtpt[1] / jtpt[0] : 0.);
    h->hr31->Fill(has3 ? jtpt[2] / jtpt[0] : 0.);
    h->hr32->Fill(has32 ? jtpt[2] / jtpt[1] : 0.);
    if (has2) h->pr21->Fill(npvgood, has2 ? jtpt[1] / jtpt[0] : 0.);
    if (has3) h->pr31->Fill(npvgood, has3 ? jtpt[2] / jtpt[0] : 0.);
    if (has32) h->pr32->Fill(npvgood, has3 ? jtpt[2] / jtpt[1] : 0.);
    h->px21->Fill(npvgood, has2 ? 1 : 0);
    h->px31->Fill(npvgood, has3 ? 1 : 0);
    h->px32->Fill(npvgood, has32 ? 1 : 0);
  }


  if (_debug) cout << "Entering jet loop" << endl << flush;


  for (int i = 0; i != njt && _pass; ++i) {

    if (_debug) cout << "Loop over jet " << i << "/" << njt << endl << flush;

    // adapt variable names from different trees
    double pt = jtpt[i];
    double eta = jteta[i];
    double energy = jte[i];
    double mass = sqrt(fabs(pow(energy,2) - pow(pt*cosh(eta),2)));
    double y = jty[i];
    double phi = jtphi[i];
    double jec = jtjesnew[i];
    bool id = _jetids[i];

    double jec2 = jtjesnew[i]/jtjes[i];

    // Tag-and-probe for composition:
    // tag in barrel and fires trigger, probe in eta bin unbiased
    // only two leading jets back-to-back, third has less than 0.3*tag pT
    if (i<2 && njt>=2 && pt>_recopt && fabs(y) >= h->ymin && fabs(y) < h->ymax){
      
      int iref = (i==0 ? 1 : 0);
      double yref = jty[iref];
      double ptref = jtpt[iref];
      //double dphi = delta_phi(y-yref, phi-jtphi[iref]); // BUG!!
      double dphi = delta_phi(phi, jtphi[iref]);
      double pt3 = (njt>=3 ? jtpt[2] : 0.);
      //if (fabs(yref) < 1.3 && dphi > 2.7 && pt3 < 0.3*ptref) { // LOOSE!

      // Find Calo match for ptref so we can match to trigger
      int icalo = -1;
      double drmin = 999.;
      for (int j = 0; j != c_njt; ++j) {
	double dr = tools::oplus(delta_phi(c_jtphi[j], jtphi[iref]),
				 fabs(c_jteta[j] - jteta[iref]));
	if (dr < drmin) {
	  icalo = j;
	  drmin = dr;
	}
      }
      double ptrefcalo = (drmin<0.5 ? c_jtpt[icalo] : 0.);

      // Find Calo match for ptprobe so we can determine trigger efficiency
      icalo = -1;
      drmin = 999.;
      for (int j = 0; j != c_njt; ++j) {
	double dr = tools::oplus(delta_phi(c_jtphi[j], jtphi[i]),
				 fabs(c_jteta[j] - jteta[i]));
	if (dr < drmin) {
	  icalo = j;
	  drmin = dr;
	}
      }
      double ptprobecalo = (drmin<1.0 ? c_jtpt[icalo] : 0.);
      double emfprobe = (drmin<1.0 ? c_jtemf[icalo] : 0.);

      if (evtid && id && _jetids[iref] &&
	  fabs(yref) < 1.3 && dphi > 2.7 && pt3 < 0.3*ptref &&
	  ptrefcalo>h->pttrg) { // trigger match
	  //(ptrefcalo>30. || ptref>56.)) { // trigger match
	
        assert(h->ptrigefftp);
	h->ptrigefftp->Fill(jtpt[i], ptprobecalo>h->pttrg ? 1. : 0.);
	h->ptrigefftp0->Fill(jtpt[i], ptprobecalo>h->pttrg ? 1. : 0.);
	// Try also using ptref as reference pt
	h->ptrigefftpref->Fill(ptref, ptprobecalo>h->pttrg ? 1. : 0.);
	h->ptrigefftpref0->Fill(ptref, ptprobecalo>h->pttrg ? 1. : 0.);

	assert(h->pncandtp);
	h->pncandtp->Fill(ptref, jtn[i], _w);
	assert(h->pnchtp);
	h->pnchtp->Fill(ptref, jtnch[i], _w);
	assert(h->pnnetp);
	h->pnnetp->Fill(ptref, jtnne[i], _w);
	assert(h->pnnhtp);
	h->pnnhtp->Fill(ptref, jtnnh[i], _w);
	assert(h->pncetp);
	h->pncetp->Fill(ptref, jtnce[i], _w);
	assert(h->pnmutp);
	h->pnmutp->Fill(ptref, jtnmu[i], _w);
	//
	assert(h->pchftp);
	h->pchftp->Fill(ptref, jtchf[i], _w);
	assert(h->pneftp);
	h->pneftp->Fill(ptref, jtnef[i], _w);
	assert(h->pnhftp);
	h->pnhftp->Fill(ptref, jtnhf[i], _w);
	assert(h->pceftp);
	h->pceftp->Fill(ptref, jtcef[i], _w);
	assert(h->pmuftp);
	h->pmuftp->Fill(ptref, jtmuf[i], _w);
	assert(h->pbetatp);
	h->pbetatp->Fill(ptref, jtbeta[i], _w);
	assert(h->pbetastartp);
	h->pbetastartp->Fill(ptref, jtbetastar[i], _w);
	//
	assert(h->pchftp2);
	h->pchftp2->Fill(ptref, jtchf2[i], _w);
	assert(h->pneftp2);
	h->pneftp2->Fill(ptref, jtnef2[i], _w);
	assert(h->pnhftp2);
	h->pnhftp2->Fill(ptref, jtnhf2[i], _w);
	assert(h->pceftp2);
	h->pceftp2->Fill(ptref, jtcef2[i], _w);
	assert(h->pmuftp2);
	h->pmuftp2->Fill(ptref, jtmuf2[i], _w);
	//
	if (drmin<0.5) h->pemftp->Fill(ptref, emfprobe);
	if (drmin<0.5) h->pcalopftp->Fill(ptref, ptprobecalo/jtpt[i]);
	if (drmin<0.5) h->pcalotp->Fill(ptref, ptprobecalo/ptref);
	if (drmin<0.5) h->ppftp->Fill(ptref, jtpt[i]/ptref);
	if (drmin<0.25) h->pemftp25->Fill(ptref, emfprobe);
	if (drmin<0.25) h->pcalopftp25->Fill(ptref, ptprobecalo/jtpt[i]);
	if (drmin<0.25) h->pcalotp25->Fill(ptref, ptprobecalo/ptref);
	if (drmin<0.25) h->ppftp25->Fill(ptref, jtpt[i]/ptref);

	if (ptref > h->ptmin && ptref < h->ptmax) {

	  h->hncandtp->Fill(jtn[i], _w);
	  h->hnchtp->Fill(jtnch[i], _w);
	  h->hnnetp->Fill(jtnne[i], _w);
	  h->hnnhtp->Fill(jtnnh[i], _w);
	  h->hncetp->Fill(jtnce[i], _w);
	  h->hnmutp->Fill(jtnmu[i], _w);
	  h->hchftp->Fill(jtchf[i], _w);
	  h->hneftp->Fill(jtnef[i], _w);
	  h->hnhftp->Fill(jtnhf[i], _w);
	  h->hceftp->Fill(jtcef[i], _w); 
	  h->hmuftp->Fill(jtmuf[i], _w); 
	  h->hbetatp->Fill(jtbeta[i], _w); 
	  h->hbetastartp->Fill(jtbetastar[i], _w); 
	  //
	  if (drmin<0.5) h->hcalopftp->Fill(ptprobecalo/jtpt[i]);
	  if (drmin<0.5) h->hcalotp->Fill(ptprobecalo/ptref);
	  if (drmin<0.5) h->hpftp->Fill(jtpt[i]/ptref);
	  //
	  assert(h->pncandtp_vsnpv);
	  h->pncandtp_vsnpv->Fill(npvgood, jtn[i], _w);
	  assert(h->pnchtp_vsnpv);
	  h->pnchtp_vsnpv->Fill(npvgood, jtnch[i], _w);
	  assert(h->pnnetp_vsnpv);
	  h->pnnetp_vsnpv->Fill(npvgood, jtnne[i], _w);
	  assert(h->pnnhtp_vsnpv);
	  h->pnnhtp_vsnpv->Fill(npvgood, jtnnh[i], _w);
	  assert(h->pncetp_vsnpv);
	  h->pncetp_vsnpv->Fill(npvgood, jtnce[i], _w);
	  assert(h->pnmutp_vsnpv);
	  h->pnmutp_vsnpv->Fill(npvgood, jtnmu[i], _w);
	  //
	  assert(h->pchftp_vsnpv);
	  h->pchftp_vsnpv->Fill(npvgood, jtchf[i], _w);
	  assert(h->pneftp_vsnpv);
	  h->pneftp_vsnpv->Fill(npvgood, jtnef[i], _w);
	  assert(h->pnhftp_vsnpv);
	  h->pnhftp_vsnpv->Fill(npvgood, jtnhf[i], _w);
	  assert(h->pceftp_vsnpv);
	  h->pceftp_vsnpv->Fill(npvgood, jtcef[i], _w);
	  assert(h->pmuftp_vsnpv);
	  h->pmuftp_vsnpv->Fill(npvgood, jtmuf[i], _w);
	  assert(h->pbetatp_vsnpv);
	  h->pbetatp_vsnpv->Fill(npvgood, jtbeta[i], _w);
	  assert(h->pbetastartp_vsnpv);
	  h->pbetastartp_vsnpv->Fill(npvgood, jtbetastar[i], _w);
	  //
	  assert(h->pchftp2_vsnpv);
	  h->pchftp2_vsnpv->Fill(npvgood, jtchf2[i], _w);
	  assert(h->pneftp2_vsnpv);
	  h->pneftp2_vsnpv->Fill(npvgood, jtnef2[i], _w);
	  assert(h->pnhftp2_vsnpv);
	  h->pnhftp2_vsnpv->Fill(npvgood, jtnhf2[i], _w);
	  assert(h->pceftp2_vsnpv);
	  h->pceftp2_vsnpv->Fill(npvgood, jtcef2[i], _w);
	  assert(h->pmuftp2_vsnpv);
	  h->pmuftp2_vsnpv->Fill(npvgood, jtmuf2[i], _w);
	  //
	  if (drmin<0.5) h->pemftp_vsnpv->Fill(npvgood, emfprobe);
	  if (drmin<0.5) h->pcalopftp_vsnpv->Fill(npvgood, ptprobecalo/jtpt[i]);
	}
      } // dijet system
    } // tag-and-probe

    // Check effect of ID cuts
    if (pt>_recopt && fabs(y) >= h->ymin && fabs(y) < h->ymax) {

      if (_debug) {
	cout << "..." << h->trigname << " | " << " index " << i << "/" << njt
	     << " jet pt: " << pt << " y : " << y
	     << " id " << id << " jec: " << jec << endl;
	cout << "...evt id: " << evtid << " weight: " << _w
	     << " met: " << met << " metsumet: " << metsumet << endl;
      }
    
      assert(h->hpt_noid);
      h->hpt_noid->Fill(pt, _w);
      assert(h->hpt_nojetid);
      if (evtid) h->hpt_nojetid->Fill(pt, _w);
      assert(h->hpt_noevtid);
      if (id)    h->hpt_noevtid->Fill(pt, _w);
      // Same versus generator pT as MC extra
      // to decouple efficiency from JEC and JER
      if (h->ismc) {
	h->hpt_noid_g->Fill(jtgenpt[i], _w);
	if (evtid) h->hpt_nojetid_g->Fill(jtgenpt[i], _w);
	if (id)    h->hpt_noevtid_g->Fill(jtgenpt[i], _w);
      }
    } // ID cuts

    // Check effect of reco y vs gen y binning
    if (h->ismc) {
      double ygen = jtgeny[i]; // use jtgeny, if available
      // GenJets matched to good reco jets in good events
      if (evtid && id && pt>_recopt && jtgenr[i] < 0.25 &&
	  fabs(ygen) >= h->ymin && fabs(ygen) < h->ymax) {
	h->hpt_gg->Fill(jtgenpt[i], _w);
      }
      // GenJets matched to any reco jets in any events
      if (pt>_recopt && jtgenr[i] < 0.25 &&
	  fabs(ygen) >= h->ymin && fabs(ygen) < h->ymax) {
	h->hpt_gg0->Fill(jtgenpt[i], _w);
      }
    }

    // Debugging JEC
    /*
    if (h->ismc && h->ymin==0 && h->ymax==0.5) {

      int igen = -1;
      // Find matching genjet
      double drmin = 999;
      for (int j = 0; j != gen_njt; ++j) {
	double dr = tools::oplus(delta_phi(gen_jtphi[j], jtphi[i]),
				 fabs(gen_jteta[j] - jteta[i]));
	if (dr < drmin) {
	  igen = j;
	  drmin = dr;
	}
      }

      if (npvgood==1 && fabs(jteta[i])<1.3 &&
	  jtgenpt[i]>=80 && jtgenpt[i]<120.) {
	cout << "jet " << i << " eta = " << jteta[i]
	     << " y " << jty[i] << " ygen = " << jtgeny[i]
	     << " pt = " << jtpt[i]
	     << " ptuncorr = " << jtptu[i]
	     << " (ratio = " << jtpt[i]/jtptu[i]
	     << " , jec = " << jtjesnew[i]
	     << " , old = " << jtjes[i] << ")"
	     << " npvgood = " << npvgood
	     << " jtgenpt = " << jtgenpt[i]
	     << " dr = " << jtgenr[i]
	     << " genpt = " << (igen!=-1 ? gen_jtpt[igen] : -1)
	     << " drmin = " << drmin
	     << " gen/reco = " << jtgenpt[i]/jtpt[i] << endl; 
      }
    } // debugging JEC
    */

    // calculate efficiencies and fill histograms
    if (evtid && id && pt>_recopt && fabs(y) >= h->ymin && fabs(y) < h->ymax) {

      if (_debug) cout << "..jec uncertainty" << endl << flush;
	
      // Get JEC uncertainty
      double unc = 0.01; // default for MC
      if (_jecUnc) {
	_jecUnc->setJetEta(eta);
	_jecUnc->setJetPt(pt);
	unc = _jecUnc->getUncertainty(true);
	//_jecUnc2->Rjet(pt, unc); // use Fall10 absolute scale uncertainty
      }

      // retrieve event-wide variables
      double dphi = (njt>1 ? delta_phi(jtphi[0], jtphi[1]) : 0.);
      double dpt = (njt>1 ? fabs(jtpt[0]-jtpt[1])/(jtpt[0]+jtpt[1]) : 0.999);
      //double met = this->met;
      //double metphi = this->metphi;
      double sumet = this->metsumet;

      // calculate and/or retrieve efficiencies
      double ideff = 1.;
      double vtxeff = 1.;
      double dqmeff = 1.;
      double trigeff = 1.;
      double eff = ideff * vtxeff * dqmeff * trigeff;
    
      if (_debug) cout << "..raw spectrum" << endl << flush;

      // For trigger efficiency
      if (h->ismc) {
	h->hpt_jt30->Fill(pt, _w0 * _wt["jt30"]);
	h->hpt_jt60->Fill(pt, _w0 * _wt["jt60"]);
	//h->hpt_jt80->Fill(pt, _w0 * _wt["jt80"]);
	h->hpt_jt110->Fill(pt, _w0 * _wt["jt110"]);
	//h->hpt_jt150->Fill(pt, _w0 * _wt["jt150"]);
	h->hpt_jt190->Fill(pt, _w0 * _wt["jt190"]);
	h->hpt_jt240->Fill(pt, _w0 * _wt["jt240"]);
	h->hpt_jt300->Fill(pt, _w0 * _wt["jt300"]);
	h->hpt_jt370->Fill(pt, _w0 * _wt["jt370"]);
      }

      // raw spectrum
      assert(h->hpt);
      h->hpt->Fill(pt, _w);
      assert(h->hpt0);
      h->hpt0->Fill(pt, _w);
      //if (y>=0) {
	//h->hpt_plus->Fill(pt, _w);
	//h->hpt0_plus->Fill(pt, _w);
	//h->hpt_plus_38x->Fill(pt, _w);
      //}
      //if (y<0) {
	//h->hpt_minus->Fill(pt, _w);
	//h->hpt0_minus->Fill(pt, _w);
	//h->hpt_minus_38x->Fill(pt, _w);
      //}
      // Do proper event statistics
      if (h->hpttmp->GetBinContent(h->hpttmp->FindBin(pt))==0)
	h->hptevt->Fill(pt, _w);
      h->hpttmp->Fill(pt);

      // leading and non-leading jets
      assert(h->hpt1);
      if (i==0) { h->hpt1->Fill(pt, _w); }
      assert(h->hpt2);
      if (i==1) { h->hpt2->Fill(pt, _w); }
      assert(h->hpt3);
      if (i==2) { h->hpt3->Fill(pt, _w); }

      if (_debug) cout << "..basic properties" << endl << flush;

      // basic properties
      h->ppt->Fill(pt, pt, _w); assert(h->ppt);
      h->pmass->Fill(pt, mass/energy, _w); assert(h->pmass);
      h->pjec->Fill(pt, jec, _w); assert(h->pjec);
      h->pjec2->Fill(pt, jec2, _w); assert(h->pjec2);
      h->punc->Fill(pt, unc, _w); assert(h->punc);

      // JEC monitoring
      h->pjec_l1->Fill(pt, jtjes_l1[i], _w); assert(h->pjec_l1);
      h->pjec_l2l3->Fill(pt, jtjes_l2l3[i], _w); assert(h->pjec_l2l3);
      h->pjec_res->Fill(pt, jtjes_res[i], _w); assert(h->pjec_res);
      
      // Pile-up information
      h->pa->Fill(pt, jta[i], _w);
      h->prho->Fill(pt, rho, _w);
      h->pnpv->Fill(pt, npvgood, _w);
      h->pnpvall->Fill(pt, npv, _w);
      if (npv==1) h->prho0->Fill(pt, rho, _w);
      if (npvgood==1) h->prho1->Fill(pt, rho, _w);
      if (npvgood==2) h->prho2->Fill(pt, rho, _w);
      if (npvgood==3) h->prho3->Fill(pt, rho, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->prhovsnpv->Fill(npvgood, rho, _w);
	h->prhovsnpvall->Fill(npv, rho, _w);
	h->h2rhovsnpv->Fill(npvgood, rho, _w);
      }
      //
      double ptraw = jtptu[i];
      double eraw = jteu[i];
      double chs = ptraw*jtchf[i]*jtbetastar[i];
      double chs_x = ptraw*jtchf[i]*(1-jtbeta[i]);
      h->pchs->Fill(pt, chs, _w);
      h->pchsx->Fill(pt, chs_x, _w);
      if (npv==1) h->pchs0->Fill(pt, chs, _w);
      if (npvgood==1) h->pchs1->Fill(pt, chs, _w);
      if (npvgood==2) h->pchs2->Fill(pt, chs, _w);
      if (npvgood==3) h->pchs3->Fill(pt, chs, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->pchsxvsnpv->Fill(npvgood, chs_x, _w);
	h->pchsvsnpv->Fill(npvgood, chs, _w);
	h->pchsvsnpvall->Fill(npv, chs, _w);
	h->h2chsvsnpv->Fill(npvgood, chs, _w);
      }
      //
      _L1Off->setJetEta(eta);
      //_L1Off->setJetEta(max(-4.7, min(4.7, eta)));
      _L1Off->setJetE(eraw);//ptraw*cosh(eta));//energy);
      _L1Off->setJetPt(ptraw);
      _L1Off->setNPV(npvgood);
      double l1off = (1-_L1Off->getCorrection())*ptraw;
      double l2l3 = jtjesnew[i] / l1off;
      h->hpt_l1off->Fill(l2l3*(ptraw-l1off), _w);
      h->pl1off->Fill(pt, l1off, _w);
      if (npv==1) h->pl1off0->Fill(pt, l1off, _w);
      if (npvgood==1) h->pl1off1->Fill(pt, l1off, _w);
      if (npvgood==2) h->pl1off2->Fill(pt, l1off, _w);
      if (npvgood==3) h->pl1off3->Fill(pt, l1off, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->pl1offvsnpv->Fill(npvgood, l1off, _w);
	h->pl1offvsnpvall->Fill(npv, l1off, _w);
      }
      //
      _L1Fast->setJetEta(eta);
      //_L1Fast->setJetEta(max(-4.7, min(4.7, eta)));
      _L1Fast->setJetPt(ptraw);
      _L1Fast->setJetE(eraw);
      _L1Fast->setJetA(jta[i]);
      _L1Fast->setRho(rho);
      double l1 = _L1Fast->getCorrection();
      //double l2l3 = jtjesnew[i] / l1;
      double l1fast = (1 - l1)*ptraw;
      h->hpt_l1fast->Fill(l2l3*(ptraw-l1fast), _w);
      h->pl1fast->Fill(pt, l1fast, _w);
      if (npv==1) h->pl1fast0->Fill(pt, l1fast, _w);
      if (npvgood==1) h->pl1fast1->Fill(pt, l1fast, _w);
      if (npvgood==2) h->pl1fast2->Fill(pt, l1fast, _w);
      if (npvgood==3) h->pl1fast3->Fill(pt, l1fast, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->pl1fastvsnpv->Fill(npvgood, l1fast, _w);
	h->pl1fastvsnpvall->Fill(npv, l1fast, _w);
      }

      // Same for CHS
      _L1OffCHS->setJetEta(eta);
      //_L1OffCHS->setJetEta(max(-4.7, min(4.7, eta)));
      _L1OffCHS->setJetE(eraw);//ptraw*cosh(eta));//energy);
      _L1OffCHS->setJetPt(ptraw);
      _L1OffCHS->setNPV(npvgood);
      double l1offchs = (1-_L1OffCHS->getCorrection())*ptraw;
      //double l2l3 = jtjesnew[i] / l1offchs;
      //h->hpt_l1offchs->Fill(l2l3*(ptraw-l1offchs), _w);
      //h->pl1offchs->Fill(pt, l1offchs, _w);
      //if (npv==1) h->pl1offchs0->Fill(pt, l1offchs, _w);
      //if (npvgood==1) h->pl1offchs1->Fill(pt, l1offchs, _w);
      //if (npvgood==2) h->pl1offchs2->Fill(pt, l1offchs, _w);
      //if (npvgood==3) h->pl1offchs3->Fill(pt, l1offchs, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->pl1offchsvsnpv->Fill(npvgood, l1offchs, _w);
	h->pl1offchsvsnpvall->Fill(npv, l1offchs, _w);
      }
      //
      _L1FastCHS->setJetEta(eta);
      //_L1FastCHS->setJetEta(max(-4.7, min(4.7, eta)));
      _L1FastCHS->setJetPt(ptraw);
      _L1FastCHS->setJetE(eraw);
      _L1FastCHS->setJetA(jta[i]);
      _L1FastCHS->setRho(rho);
      double l1chs = _L1FastCHS->getCorrection();
      //double l2l3 = jtjesnew[i] / l1;
      double l1fastchs = (1 - l1chs)*ptraw;
      //h->hpt_l1fastchs->Fill(l2l3*(ptraw-l1fastchs), _w);
      //h->pl1fastchs->Fill(pt, l1fastchs, _w);
      //if (npv==1) h->pl1fastchs0->Fill(pt, l1fastchs, _w);
      //if (npvgood==1) h->pl1fastchs1->Fill(pt, l1fastchs, _w);
      //if (npvgood==2) h->pl1fastchs2->Fill(pt, l1fastchs, _w);
      //if (npvgood==3) h->pl1fastchs3->Fill(pt, l1fastchs, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
	h->pl1fastchsvsnpv->Fill(npvgood, l1fastchs, _w);
	h->pl1fastchsvsnpvall->Fill(npv, l1fastchs, _w);
      }


      // efficiencies
      assert(h->peff);
      h->peff->Fill(pt, eff, _w);
      assert(h->pideff);
      h->pideff->Fill(pt, ideff, _w);
      assert(h->pvtxeff);
      h->pvtxeff->Fill(pt, vtxeff, _w);
      assert(h->pdqmeff);
      h->pdqmeff->Fill(pt, dqmeff, _w);
      assert(h->ptrigeff);
      h->ptrigeff->Fill(pt, trigeff, _w);

      if (_debug) cout << "..control plots of components" << endl << flush;

      // control plots of jet components (JEC)
      assert(h->pncand);
      h->pncand->Fill(pt, jtn[i], _w);
      assert(h->pnch);
      h->pnch->Fill(pt, jtnch[i], _w);
      assert(h->pnne);
      h->pnne->Fill(pt, jtnne[i], _w);
      assert(h->pnnh);
      h->pnnh->Fill(pt, jtnnh[i], _w);
      assert(h->pnce);
      h->pnce->Fill(pt, jtnce[i], _w);
      assert(h->pnmu);
      h->pnmu->Fill(pt, jtnmu[i], _w);
      //
      assert(h->pchf);
      h->pchf->Fill(pt, jtchf[i], _w);
      assert(h->pnef);
      h->pnef->Fill(pt, jtnef[i], _w);
      assert(h->pnhf);
      h->pnhf->Fill(pt, jtnhf[i], _w);
      assert(h->pcef);
      h->pcef->Fill(pt, jtcef[i], _w);
      assert(h->pmuf);
      h->pmuf->Fill(pt, jtmuf[i], _w);
      assert(h->pbeta);
      h->pbeta->Fill(pt, jtbeta[i], _w);
      assert(h->pbetastar);
      h->pbetastar->Fill(pt, jtbetastar[i], _w);
      
      // Find Calo match for ptprobe so we can determine trigger efficiency
      int icalo = -1;
      double drmin = 999.;
      for (int j = 0; j != c_njt; ++j) {
	double dr = tools::oplus(delta_phi(c_jtphi[j], phi),
				 fabs(c_jteta[j] - eta));
	if (dr < drmin) {
	  icalo = j;
	  drmin = dr;
	}
      }
      double ptcalo = (drmin<0.5 ? c_jtpt[icalo] : 0.);
      double emf = (drmin<0.5 ? c_jtemf[icalo] : 0.);
      if (drmin<0.5) h->pemf->Fill(pt, emf);
      if (drmin<0.25) h->pemf25->Fill(pt, emf);
      if (drmin<0.5) h->pcalopf->Fill(pt, ptcalo/pt);
      if (drmin<0.25) h->pcalopf25->Fill(pt, ptcalo/pt);

      // for b-jets

      // control plots for topology (JEC)
      if (pt >= h->ptmin && pt < h->ptmax) {

	if (_debug) cout << "..control plots for topology" << endl << flush;

	if (h->ismc) {
	  h->hitpu->Fill(itpu, _w);
	  h->hootpuearly->Fill(ootpuearly, _w);
	  h->hootpulate->Fill(ootpulate, _w);
	  h->h2itvsoot->Fill(itpu, ootpulate, _w);
	}

	h->hnpvgood->Fill(npvgood, _w);
	h->hselpt->Fill(pt, _w);
	h->hmass->Fill(mass/energy, _w);
	h->hy->Fill(y, _w);
	h->hy2->Fill(y, _w);
	h->heta->Fill(eta, _w);
	h->heta2->Fill(eta, _w);
	h->hphi->Fill(phi, _w);
	h->hdphi->Fill(dphi, _w);
	h->hdpt->Fill(dpt, _w);
	h->hjet->Fill(pt / sumet, _w);
	h->hmet->Fill(met / sumet, _w);
	h->hmetphi->Fill(delta_phi(metphi, phi), _w);
	// control plots for vertex
	h->hpvndof->Fill(pvndof);
	h->hpvx->Fill(pvx-bsx);
	h->hpvy->Fill(pvy-bsy);
	h->hpvz->Fill(pvz-0.);
	h->hpvr->Fill(tools::oplus(pvx-bsx, pvy-bsy));
	h->hpvrho->Fill(pvrho-tools::oplus(bsx, bsy));
	// closure plots for JEC
	h->hmpf->Fill(1 + met * cos(delta_phi(metphi, phi)) / pt, _w);
	h->hmpf1->Fill(1 + met1 * cos(delta_phi(metphi1, phi)) / pt, _w);
	h->hmpf2->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);
	//
	if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2] < _recopt) &&
	    (njt>=2 && jtpt[1] > _recopt && delta_phi(jtphi[0],jtphi[1]) > 2.7))
	  h->hmpfy->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);

	// Component fractions
	h->hncand->Fill(jtn[i], _w);
	h->hnch->Fill(jtnch[i], _w);
	h->hnne->Fill(jtnne[i], _w);
	h->hnnh->Fill(jtnnh[i], _w);
	h->hnce->Fill(jtnce[i], _w);
	h->hnmu->Fill(jtnmu[i], _w);
	//
	h->hchf->Fill(jtchf[i], _w);
	h->hnef->Fill(jtnef[i], _w);
	h->hnhf->Fill(jtnhf[i], _w);
	h->hcef->Fill(jtcef[i], _w);
	h->hmuf->Fill(jtmuf[i], _w);
	h->hbeta->Fill(jtbeta[i], _w);
	h->hbetastar->Fill(jtbetastar[i], _w);
	//
	if (drmin<0.5) h->hcalopf->Fill(ptcalo/pt); // fix nov28

	h->hyeta->Fill(TMath::Sign(y-eta,y), _w);
	h->hyeta2->Fill(y-eta, _w);
	h->hbetabetastar->Fill(jtbeta[i], jtbetastar[i], _w);
	h->hetaphi->Fill(eta, phi, _w);
	if (_rd) {
	  if (run<=167913)
	    h->hetaphiA->Fill(eta, phi, _w);
	  else
	    h->hetaphiB->Fill(eta, phi, _w);
	}
	
      } // within trigger pT range
      
      int iprobe = i;
      int itag = (iprobe==0 ? 1 : 0);
      double pttag = (njt>=2 ? jtpt[itag] : 0);
      if (iprobe<2 && pttag >= h->ptmin && pttag < h->ptmax &&
	  fabs(jty[itag]) < 1.3) {

	if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2] < _recopt) &&
	    (njt>=2 && jtpt[1] > _recopt && delta_phi(jtphi[0],jtphi[1]) > 2.7))
	  h->hmpfx->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pttag, _w);
      }
      double ptave = (njt>=2 ? 0.5*(pt+pttag) : 0);
      if (iprobe<2 && ptave >= h->ptmin && ptave < h->ptmax &&
	  fabs(jty[itag]) < 1.3) {

	if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2] < _recopt) &&
	    (njt>=2 && jtpt[1] > _recopt && delta_phi(jtphi[0],jtphi[1]) > 2.7))
	  h->hmpfz->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / ptave, _w);
      }

      // closure plots for JEC
      h->pdpt->Fill(pt, dpt, _w);
      h->pmpf->Fill(pt, 1 + met * cos(delta_phi(metphi, phi)) / pt, _w);
      h->pmpf1->Fill(pt, 1 + met1 * cos(delta_phi(metphi1, phi)) / pt, _w);
      h->pmpf2->Fill(pt, 1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);
      //
      if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2] < _recopt) &&
	  (njt>=2 && jtpt[1] > _recopt && delta_phi(jtphi[0],jtphi[1]) > 2.7)) {

	h->pmpfx->Fill(pttag, 1 + met2 * cos(delta_phi(metphi2, phi))/pttag, _w);
	h->pmpfy->Fill(pt, 1 + met2 * cos(delta_phi(metphi2, phi))/pt, _w);
	h->pmpfz->Fill(ptave, 1 + met2 * cos(delta_phi(metphi2, phi))/ptave, _w);
      }

      // MC extras
      if (_mc && jtgenr[i]<0.25) {
	h->hpt_gtw->Fill(jtgenpt[i], _w);
      }
      //
      if (h->ismc) {

	if (jtgenr[i]<0.25) {

	  //int flv = abs(jtgenflv[i]);
	  double ptgen = jtgenpt[i];
	  //double r = (jtgenpt[i] ? pt/jtgenpt[i] : 0);
	  double r = (ptgen ? pt/ptgen : 0);
	  //double resp = (jtjesnew[i] ? r / jtjesnew[i] : 0);
	  double dy = (r ? TMath::Sign(jty[i]-jtgeny[i], jtgeny[i]) : 0.);
	  h->hpt_r->Fill(pt, _w);
	  h->hpt_g->Fill(ptgen, _w);
	  h->ppt_r->Fill(pt, pt, _w);
	  h->ppt_g->Fill(ptgen, ptgen, _w);

	  // Response closure vs NPV
	  //if (r) h->p3rvsnpv->Fill(ptgen, jteta[i], npvgood, resp, _w);
	  //if (r) h->p3rvsnpvW->Fill(ptgen, fabs(jteta[i]), npvgood, resp, _w);
	  if (r) h->p2rvsnpv->Fill(ptgen, npvgood, r, _w);

	  // Response closure
	  if (r) h->h2r_r->Fill(pt, r, _w);
	  if (r) h->h2r_g->Fill(ptgen, r, _w);
	  if (r) h->p2r_r->Fill(pt, r, _w);
	  if (r) h->p2r_g->Fill(ptgen, r, _w);
	  if (r) h->p2r_ruw->Fill(pt, r); // unweighted!
	  if (r) h->p2r_guw->Fill(ptgen, r); // unweighted!
	  
	  // Rapidity closure
	  if (r) h->h2dy_r->Fill(pt, dy, _w);
	  if (r) h->h2dy_g->Fill(ptgen, dy, _w); 
	  if (r) h->p2dy_r->Fill(pt, dy, _w);
	  if (r) h->p2dy_g->Fill(ptgen, dy, _w); 
	  if (r) h->p2dy_ruw->Fill(pt, dy); // unweighted
	  if (r) h->p2dy_guw->Fill(ptgen, dy); // unweighted
	  if (r) h->pdy_r->Fill(pt, fabs(y), dy, _w);
	  if (r) h->pdy_g->Fill(ptgen, fabs(y), dy, _w);
	}
      } // is MC

    } // if id && etabin

    // Fill outside of eta bin also
    if (h->ismc && evtid && id && pt > _recopt &&
	jtgenr[i]<0.25 && jtgenpt[i]!=0 && jtjesnew[i]!=0) {

      double ptgen = jtgenpt[i];
      double r = pt / ptgen;
      double resp = r / jtjesnew[i];

      // Response closure vs NPV
      h->p3rvsnpv->Fill(ptgen, jteta[i], npvgood, resp, _w);
      h->p3rvsnpvW->Fill(ptgen, fabs(jteta[i]), npvgood, resp, _w);
    } // if id && MC
  } // for i

  if (h->dofriends) {
    
    bool caloevtid = (c_met < 0.5 * c_metsumet || c_met < 50.); // 100->50
    
    for (int i = 0; i != c_njt; ++i) {
      
      double y = c_jty[i];
      if (caloevtid && c_jtidtight[i] && _pass &&
	  fabs(y) >= h->ymin && fabs(y) < h->ymax) {
	
        double pt = c_jtpt[i];
	
        if (pt > _recopt) {
	  
          h->hpt_ak5calo->Fill(pt, _w);

	  // JEC monitoring
	  h->pcjec_l1->Fill(pt, c_jtjes_l1[i], _w); assert(h->pcjec_l1);
	  h->pcjec_l2l3->Fill(pt, c_jtjes_l2l3[i], _w); assert(h->pcjec_l2l3);
	  h->pcjec_res->Fill(pt, c_jtjes_res[i], _w); assert(h->pcjec_res);

	  //h->pa_ak5calo->Fill(pt, c_jta[i], _w);
	  h->pa_ak5calo->Fill(pt, c_jta, _w);
          h->prho_ak5calo->Fill(pt, c_rho, _w);
          if (npv==1) h->prho0_ak5calo->Fill(pt, c_rho, _w);
          if (npvgood==1) h->prho1_ak5calo->Fill(pt, c_rho, _w);
          if (npvgood==2) h->prho2_ak5calo->Fill(pt, c_rho, _w);
          if (npvgood==3) h->prho3_ak5calo->Fill(pt, c_rho, _w);
          if (pt >= h->ptmin && pt < h->ptmax) {
            h->prhovsnpv_ak5calo->Fill(npvgood, c_rho, _w);
            h->prhovsnpvall_ak5calo->Fill(npv, c_rho, _w);
            h->h2rhovsnpv_ak5calo->Fill(npvgood, c_rho, _w);
          }
          //                                                                    
          double ptraw = c_jtptu[i];
          double eraw = c_jtptu[i];
	  //
          _L1Off_ak5calo->setJetEta(c_jteta[i]);//y);
	  //_L1Off_ak5calo->setJetEta(max(-4.7, min(4.7, y)));
          _L1Off_ak5calo->setJetPt(ptraw);
          _L1Off_ak5calo->setJetE(eraw);//ptraw*cosh(y));
          _L1Off_ak5calo->setNPV(npvgood);
          double l1off = (1-_L1Off_ak5calo->getCorrection())*ptraw;
	  double l2l3 = c_jtjesnew[i] / l1off;
          h->hpt_l1off_ak5calo->Fill(l2l3*(ptraw-l1off), _w);
          h->pl1off_ak5calo->Fill(pt, l1off, _w);
          if (npv==1) h->pl1off0_ak5calo->Fill(pt, l1off, _w);
          if (npvgood==1) h->pl1off1_ak5calo->Fill(pt, l1off, _w);
          if (npvgood==2) h->pl1off2_ak5calo->Fill(pt, l1off, _w);
          if (npvgood==3) h->pl1off3_ak5calo->Fill(pt, l1off, _w);
          if (pt >= h->ptmin && pt < h->ptmax) {
            h->pl1offvsnpv_ak5calo->Fill(npvgood, l1off, _w);
            h->pl1offvsnpvall_ak5calo->Fill(npv, l1off, _w);
          }
          //                                                                    
          _L1Fast_ak5calo->setJetEta(c_jteta[i]);//y);
	  //_L1Fast_ak5calo->setJetEta(max(-4.7, min(4.7, y)));
          _L1Fast_ak5calo->setJetPt(ptraw);
          _L1Fast_ak5calo->setJetE(eraw);
          _L1Fast_ak5calo->setJetA(c_jta);//[i]);
          _L1Fast_ak5calo->setRho(c_rho);
          double l1 = _L1Fast_ak5calo->getCorrection();
          //double l2l3 = c_jtjesnew[i] / l1;
          double l1fast = (1 - l1)*ptraw;
          h->hpt_l1fast_ak5calo->Fill(l2l3*(ptraw-l1fast), _w);
          h->pl1fast_ak5calo->Fill(pt, l1fast, _w);
          if (npv==1) h->pl1fast0_ak5calo->Fill(pt, l1fast, _w);
          if (npvgood==1) h->pl1fast1_ak5calo->Fill(pt, l1fast, _w);
          if (npvgood==2) h->pl1fast2_ak5calo->Fill(pt, l1fast, _w);
          if (npvgood==3) h->pl1fast3_ak5calo->Fill(pt, l1fast, _w);
          if (pt >= h->ptmin && pt < h->ptmax) {
            h->pl1fastvsnpv_ak5calo->Fill(npvgood, l1fast, _w);
            h->pl1fastvsnpvall_ak5calo->Fill(npv, l1fast, _w);
          }

          if (pt >= h->ptmin && pt < h->ptmax) {
	    h->hetaphi_calo->Fill(c_jteta[i], c_jtphi[i], _w);
	    if (_rd) {
	      if (run<=167913) 
		h->hetaphiA_calo->Fill(c_jteta[i], c_jtphi[i], _w);
	      else
		h->hetaphiB_calo->Fill(c_jteta[i], c_jtphi[i], _w);
	    }
	  }
	} // reco pt
      } // y bin                                                                
    }// for i 
  } // dofriends

  // Unbiased generator spectrum
  if (_mc) {
    for (int i = 0; i != gen_njt; ++i) {
      double y = gen_jty[i];
      if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {
	h->hpt_g0tw->Fill(gen_jtpt[i], _w);
      }
    }
  }
  //
  if (h->ismc) {
    
    // unfolding studies (Mikael)
    for (int i = 0; i != gen_njt; ++i) {

      double ygen = fabs(gen_jty[i]);

      for (int j = 0; j != njt; ++j) {

	//double yreco = fabs(jty[j]);
	bool id = (_jetids[j] && evtid && _pass);
	
	if ((ygen >= h->ymin && ygen < h->ymax && gen_jtpt[i]>_recopt) &&
	    //(yreco >= h->ymin && yreco < h->ymax) 
	    (jtpt[j]>_recopt && id)) {
	
	  double dr = tools::oplus(delta_phi(gen_jtphi[i], jtphi[j]),
				   fabs(gen_jteta[i] - jteta[j]));
	  if (dr < 0.25) {
	    h->mT->Fill(gen_jtpt[i], jtpt[j], _w);
	    h->mTf->Fill(gen_jtpt[i], jtpt[j], _w);
	    h->mTuw->Fill(gen_jtpt[i], jtpt[j]);
	    h->mTfuw->Fill(gen_jtpt[i], jtpt[j]);
	  }
	} // rapidity bin
      } // for j
    } // for i
    //
    for (int i = 0; i != gen_njt; ++i) {
      double ygen = fabs(gen_jty[i]);
      if (ygen >= h->ymin && ygen < h->ymax && gen_jtpt[i]>_recopt) {
	h->mx->Fill(gen_jtpt[i], _w);
	h->mxf->Fill(gen_jtpt[i], _w);
	h->mxuw->Fill(gen_jtpt[i]);
	h->mxfuw->Fill(gen_jtpt[i]);
      }
    } // for i
    //
    for (int j = 0; j != njt; ++j) {
      double yreco = fabs(jty[j]);
      bool id  = (_jetids[j] && evtid && _pass && jtpt[j]>_recopt);
      if (yreco >= h->ymin && yreco < h->ymax && id) {
	h->my->Fill(jtpt[j], _w);
	h->myf->Fill(jtpt[j], _w);
	h->myuw->Fill(jtpt[j]);
	h->myfuw->Fill(jtpt[j]);
      }
    } // for j
    
    for (int i = 0; i != gen_njt; ++i) {
      double y = gen_jty[i];
      if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {
	h->hpt_g0->Fill(gen_jtpt[i], _w);
	//if (abs(gen_jtgenflv[i])==5) {
	//h->hbpt_g0->Fill(gen_jtpt[i], _w);
	//h->pbpt_g0->Fill(gen_jtpt[i], gen_jtpt[i], _w);
	//}
      }

      // Debugging JEC
      /*
      if (h->ymin==0 && h->ymax==0.5) {

	int ireco = -1;
	// Find matching genjet
	double drmin = 999;
	for (int j = 0; j != njt; ++j) {
	  double dr = tools::oplus(delta_phi(gen_jtphi[i], jtphi[j]),
				   fabs(gen_jteta[i] - jteta[j]));
	  if (dr < drmin) {
	    ireco = j;
	    drmin = dr;
	  }
	}
	
	if (fabs(gen_jty[i])<1.3 && gen_jtpt[i]>=80 && gen_jtpt[i]<120.) {
	  if (drmin>0.25) cout << "***** mismatch!!!" << endl;
	  cout << "*genjet " << i
	       << ", genpt = " << gen_jtpt[i]
	       << ", geny = " << gen_jty[i]
	       << ", recopt = " << (ireco!=-1 ? jtpt[ireco] : -1)
	       << ", recoy = " << (ireco!=-1 ? jty[ireco] : 0)
	       << ", drmin = " << drmin << endl;
	}
      } // debugging JEC
      */

    } // for i
  } // gen spectrum

} // fillBasic

// Write and delete histograms
void fillHistos::writeBasics() {
  
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeBasics():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (map<string, vector<basicHistos*> >::iterator it = _histos.begin();
       it != _histos.end(); ++it) {

    for (unsigned int i = 0; i != it->second.size(); ++i) {

      // Luminosity information
      basicHistos *h = it->second[i];
      for (int j = 0; j != h->hlumi->GetNbinsX()+1; ++j) {
	h->hlumi->SetBinContent(j, _rd ? h->lumsum : 1. );
	h->hlumi2->SetBinContent(j, _rd ? h->lumsum2 : 1. );
      }
      
      delete h;//it->second[i];
    } // for i
  } // for it
  
  cout << "\nOutput stored in " << _outfile->GetName() << endl;
  _outfile->Close();
  _outfile->Delete();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeBasic() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeBasic


// Initialize basic histograms for trigger and eta bins
 void fillHistos::initRunHistos(string name, double ymin, double ymax) {

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initRunHistos("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile :
	      new TFile(Form("output-%s-1.root",_type.c_str()), "RECREATE"));
  assert(f && !f->IsZombie());
  assert(f->mkdir(name.c_str()));
  assert(f->cd(name.c_str()));
  TDirectory *dir = gDirectory;

  runHistos *h = new runHistos(dir, ymin, ymax);
  _runhistos[name] = h;

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initRunHistos() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initRunHistos

// Fill run histograms
void fillHistos::fillRunHistos(string name) {

  runHistos *h = _runhistos[name];
  assert(h);

  // Luminosity information
  if (_rd && h->lums[run][lbn]==0) {

    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    if (lum==0) return; // no lumi no pass

    h->lumsum += lum;
    h->lumsum2 += lum2;
    h->runlums[run] += lum;
    h->runlums2[run] += lum2;
    h->lums[run][lbn] = 1;

    for (unsigned int i = 0; i != h->trg.size(); ++i) {

      string const& t = h->trg[i];
      double prescale(0);
      if (_prescales[t].find(run)==_prescales[t].end()) {
	if (_trigs.find(t)!=_trigs.end()) {
	  *ferr << "Prescale not found for trigger " << t
		<< " run " << run << endl << flush;
	  assert(false);
	}
      }
      else 
	prescale = _prescales[t][run]; //assert(prescale);
      h->runlums_trg[t][run] += (prescale ? lum / prescale : 0.);
    } // for i

    // Initialize counters for a new run
    if (h->lums[run].size()==1) {
      for (unsigned int i = 0; i != h->trg.size(); ++i) {

	string const& t = h->trg[i];
	h->p_trg[t][run] = 0;
	h->t_trg[t][run] = 0;
	h->npv_trg[t][run] = 0;
	h->c_chf[t][run] = 0;
	h->c_nef[t][run] = 0;
	h->c_nhf[t][run] = 0;
	h->c_betastar[t][run] = 0;
	h->t_trgtp[t][run] = 0;
	h->c_chftp[t][run] = 0;
	h->c_neftp[t][run] = 0;
	h->c_nhftp[t][run] = 0;
	h->c_betastartp[t][run] = 0;

	for (unsigned int j = 0; j != h->trg.size(); ++j) {

	  string const& t2 = h->trg[j];
	  h->p_trgpair[t+t2][run] = 0;
	} // for j
      } // for i
    } // new run
  }

  bool evtid = (met < 0.5 * metsumet || met < 100.);

  double dphi = (njt>=2 ? delta_phi(jtphi[0], jtphi[1]) : 0.);
  double pt3 = (njt>=3 ? jtpt[2] : 0.);

  for (int i = 0; i != njt; ++i) {

    double pt = jtpt[i];
    double y = jty[i];

    if (h->ymin <= fabs(y) && fabs(y) < h->ymax && _pass && _jetids[i]
	&& evtid) {// && caloevtid) {

      for (set<string>::const_iterator it = _trigs.begin(); it != _trigs.end();
	   ++it) {

	string const& t = *it;

	if (pt > 18.) ++h->p_trg[t][run];
	if (pt > h->pt[t]) {
	  ++h->t_trg[t][run];
	  h->npv_trg[t][run] += npv;
	  h->npvgood_trg[t][run] += npvgood;
	  h->c_chf[t][run] += jtchf[i];
	  h->c_nef[t][run] += jtnef[i];
	  h->c_nhf[t][run] += jtnhf[i];
	  h->c_betastar[t][run] += jtbetastar[i];
	}

	int iref = (i==0 ? 1 : 0);
	if (i<2 && dphi > 2.7 && pt3 < jtpt[iref] && fabs(jty[iref]) < 1.3 &&
	    jtpt[iref] > h->pt[t] && _jetids[iref]) {
	  ++h->t_trgtp[t][run];
	  h->c_chftp[t][run] += jtchf[i];
	  h->c_neftp[t][run] += jtnef[i];
	  h->c_nhftp[t][run] += jtnhf[i];
	  h->c_betastartp[t][run] += jtbetastar[i];
	}
	
	for (set<string>::const_iterator jt = _trigs.begin(); 
	     jt != _trigs.end(); ++jt) {
	  
	  string const& t2 = *jt;
	  if (t!=t2) ++h->p_trgpair[t+t2][run];
	} // for jt
      } // for it
    }
  } // for i

} // fillRunHistos

// Write and delete histograms
void fillHistos::writeRunHistos() {

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeRunHistos():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (map<string, runHistos*>::iterator it = _runhistos.begin();
       it != _runhistos.end(); ++it) {

    runHistos *h = it->second;
    delete h;
  } // for it
  
  cout << "\nOutput (runHistos) stored in " << _outfile->GetName() << endl;
  //_outfile->Close();
  //_outfile->Delete();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeRunHistos() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
	       info.fMemTotal, info.fMemUsed, info.fMemFree,
	       info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeRunHistos


void fillHistos::fillJetID(vector<bool> &id) {

  assert(int(id.size())==njt);

  for (int i = 0; i != njt; ++i) {

    id[i] = ((fabs(jty[i])<2.5 ? jtidtight[i] : jtidloose[i]));
	     //&& (_mc || jtbeta[i]!=0));
	     //&& (_mc || (1-jtbetastar[i])>0.5));

    if (_doECALveto) {
      assert(ecalveto);
      int ibin = ecalveto->FindBin(jteta[i],jtphi[i]);
      id[i] = (id[i] && ecalveto->GetBinContent(ibin)==0);
    }
  }

} // fillJetID


// Load good run and LS information
void fillHistos::loadJSON(const char* filename) {

  cout << "Processing loadJSON(\"" << filename << "\"..." << endl;
  ifstream file(filename, ios::in);
  assert(file.is_open());
  char c;
  string s, s2;
  char s1[256];
  int rn, ls1, ls2, nrun(0), nls(0);
  file.get(c);
  assert(c=='{');
  while (file >> s && sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (_debug)
      cout << "\"" << rn << "\": " << flush;

    while (file.get(c) && c==' ') {};
    if (_debug)
      cout << c << flush; assert(c=='[');
    ++nrun;

    bool endrun = false;
    while (!endrun && file >> s >> s2 &&
	   sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (_debug)
	cout << "["<<ls1<<","<<ls2<<"]"<<s1 << flush;

      for (int ls = ls1; ls != ls2+1; ++ls) {
	//assert(_json[rn].find(ls)==_json[rn].end()); // ok if 2 JSON files
	_json[rn][ls] = 1;
	++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (!endrun && s2!=",") {
	if (_debug)
	  cout<<"s1: "<<s2<<endl<<flush; assert(s2==",");
      }
    } // while ls
    if (_debug)
      cout << endl;

    if (s2=="]}") continue;
    else if (s2!="],") {
      if (_debug)
	cout<<"s2: "<<s2<<endl<<flush; assert(s2=="],");
    }
  } // while run
  if (s2!="]}") { cout<<"s3: "<<s2<<endl<<flush; assert(s2=="]}"); }

  cout << "Called loadJSON(\"" << filename << "\"):" << endl;
  cout << "Loaded " << nrun << " good runs and " << nls
       << " good lumi sections" << endl;

} // loadJSON


// Load luminosity information
void fillHistos::loadLumi(const char* filename) {

  cout << "Processing loadLumi(\"" << filename << "\")..." << endl;

  // Check lumi against the list of good runs
  const int a_goodruns[] = {
    // **** 4.7/fb, Nov13th ****
    // cat lumicalc/Nov13th/Runs.txt | awk '{print $1","}' | tr '\n' ' '
    // May10ReReco
    160431, 160577, 160578, 160871, 160872, 160873, 160874, 160939, 160940, 160942, 160943, 160955, 160956, 160957, 160998, 161008, 161016, 161103, 161106, 161107, 161113, 161116, 161117, 161119, 161176, 161217, 161222, 161223, 161233, 161310, 161311, 161312, 162762, 162765, 162803, 162808, 162811, 162822, 162825, 162826, 162828, 162909, 163046, 163069, 163071, 163078, 163232, 163233, 163234, 163235, 163237, 163238, 163252, 163255, 163261, 163270, 163286, 163289, 163296, 163297, 163300, 163301, 163302, 163332, 163333, 163334, 163337, 163338, 163339, 163340, 163358, 163369, 163370, 163371, 163372, 163374, 163375, 163376, 163378, 163385, 163387, 163402, 163475, 163476, 163478, 163479, 163480, 163481, 163482, 163483, 163582, 163583, 163584, 163585, 163586, 163587, 163588, 163589, 163596, 163630, 163655, 163657, 163658, 163659, 163660, 163661, 163662, 163663, 163664, 163668, 163738, 163757, 163758, 163759, 163760, 163761, 163763, 163765, 163795, 163796, 163817, 163869,
    // PromptV4
    165088, 165098, 165099, 165102, 165103, 165120, 165121, 165205, 165208, 165364, 165402, 165415, 165467, 165472, 165486, 165487, 165506, 165514, 165548, 165558, 165567, 165570, 165617, 165620, 165633, 165970, 165993, 166011, 166033, 166034, 166049, 166149, 166150, 166161, 166163, 166164, 166346, 166374, 166380, 166408, 166429, 166438, 166462, 166486, 166502, 166512, 166514, 166530, 166554, 166563, 166565, 166699, 166701, 166763, 166781, 166782, 166784, 166787, 166839, 166841, 166842, 166859, 166860, 166861, 166864, 166888, 166889, 166890, 166894, 166895, 166911, 166922, 166923, 166946, 166950, 166960, 166966, 166967, 167039, 167041, 167043, 167078, 167098, 167102, 167103, 167151, 167281, 167282, 167284, 167551, 167673, 167674, 167675, 167676, 167740, 167746, 167754, 167784, 167786, 167807, 167830, 167898, 167913,
    // Aug5ReReco    
    170826, 170842, 170854, 170876, 170896, 170899, 170901, 171050, 171091, 171098, 171102, 171106, 171117, 171156, 171178, 171219, 171274, 171282, 171315, 171369, 171446, 171484, 171578, 171812, 171876, 171880, 171895, 171897, 171921, 171926, 172014, 172024, 172033, 172163, 172208, 172252, 172254, 172255, 172268, 172286, 172389, 172399, 172400, 172401, 172411, 172478, 172619,
    // PromptV6
    172620, 172630, 172635, 172778, 172791, 172798, 172799, 172801, 172802, 172819, 172822, 172824, 172847, 172865, 172868, 172949, 172951, 172952, 172992, 172999, 173198, 173236, 173240, 173241, 173243, 173380, 173381, 173389, 173406, 173430, 173438, 173439, 173657, 173658, 173659, 173660, 173661, 173663, 173664, 173692,
    // Run2011B
    175860, 175863, 175865, 175866, 175872, 175873, 175874, 175877, 175881, 175886, 175887, 175888, 175906, 175910, 175921, 175973, 175974, 175975, 175976, 175990, 176023, 176161, 176163, 176165, 176167, 176169, 176201, 176202, 176206, 176207, 176286, 176289, 176304, 176308, 176309, 176467, 176468, 176469, 176545, 176547, 176548, 176697, 176701, 176702, 176765, 176771, 176795, 176796, 176797, 176799, 176801, 176805, 176807, 176841, 176844, 176848, 176850, 176860, 176868, 176885, 176886, 176889, 176928, 176929, 176933, 176959, 176982, 177053, 177074, 177088, 177095, 177096, 177131, 177138, 177139, 177140, 177141, 177183, 177184, 177201, 177222, 177317, 177318, 177319, 177449, 177452, 177718, 177719, 177730, 177776, 177782, 177783, 177788, 177789, 177790, 177791, 177875, 177878, 178098, 178099, 178100, 178101, 178102, 178110, 178116, 178151, 178160, 178162, 178365, 178367, 178380, 178420, 178421, 178424, 178479, 178703, 178708, 178712, 178724, 178731, 178738, 178786, 178803, 178840, 178854, 178866, 178871, 178920, 178970, 178985, 179411, 179434, 179452, 179476, 179497, 179547, 179558, 179563, 179889, 179959, 179977, 180072, 180076, 180093, 180241, 180250, 180252};
  /*
    // **** 3.2/fb, Oct11th ****
    // May10ReReco
    160431, 160577, 160578, 160871, 160872, 160873, 160874, 160939, 160940,
    160942, 160943, 160955, 160956, 160957, 160998, 161008, 161016, 161103,
    161106, 161107, 161113, 161116, 161117, 161119, 161156,
    161176, 161217, 161222,
    161223, 161233, 161310, 161311, 161312, 162762, 162765, 162803, 162808,
    162811, 162822, 162825, 162826, 162828, 162909, 163046, 163069, 163071,
    163078, 163232, 163233, 163234, 163235, 163237, 163238, 163252, 163255,
    163261, 163270, 163286, 163289, 163296, 163297, 163300, 163301, 163302,
    163332, 163333, 163334, 163337, 163338, 163339, 163340, 163358, 163369,
    163370, 163371, 163372, 163374, 163375, 163376, 163378, 163385, 163387,
    163402, 163475, 163476, 163478, 163479, 163480, 163481, 163482, 163483,
    163582, 163583, 163584, 163585, 163586, 163587, 163588, 163589, 163596,
    163630, 163655, 163657, 163658, 163659, 163660, 163661, 163662, 163663,
    163664, 163668, 163738, 163757, 163758, 163759, 163760, 163761, 163763,
    163765, 163795, 163796, 163817, 163869,
    // PromptV4
    165088, 165098, 165099, 165102,
    165103, 165120, 165121, 165205, 165208, 165364, 165402, 165415, 165467,
    165472, 165486, 165487, 165506, 165514, 165548, 165558, 165567, 165570,
    165617, 165620, 165633, 165970, 165993, 166011, 166033, 166034, 166049,
    166149, 166150, 166161, 166163, 166164, 166374, 166380, 166408, 166429,
    166438, 166462, 166486, 166502, 166512, 166514, 166530, 166554, 166563,
    166565, 166699, 166701, 166763, 166781, 166782, 166784, 166787, 166839,
    166841, 166842, 166859, 166860, 166861, 166864, 166888, 166889, 166890,
    166894, 166895, 166911, 166922, 166923, 166946, 166950, 166960, 166966,
    166967, 167039, 167041, 167043, 166346, 167078, 167098, 167102, 167103,
    167151, 167281, 167282, 167284, 167551, 167673, 167674, 167675, 167676,
    167740, 167746, 167754, 167784, 167786, 167807, 167830, 167898, 167913,
    // Aug5ReReco
    170722, 170826, 170842, 170854, 170876, 170896, 170899, 170901, 171050,
    171091, 171098, 171102, 171106, 171117, 171156, 171178, 171219, 171274,
    171282, 171315, 171369, 171446, 171484, 171578, 171812, 171876, 171880,
    171895, 171897, 171921, 171926, 172014, 172024, 172033, 172163, 172208,
    172252, 172254, 172255, 172268, 172286, 172389, 172399, 172400, 172401,
    172411, 172478, 172619, 
    // PromptV6
    172620, 172630, 172635, 172778, 172791, 172798,
    172799, 172801, 172802, 172819, 172822, 172824, 172847, 172865, 172868,
    172949, 172951, 172952, 172992, 172999, 173198, 173236, 173240, 173241,
    173243, 173380, 173381, 173389, 173406, 173430, 173438, 173439, 173657,
    173658, 173659, 173660, 173661, 173663, 173664, 173692,
    // Run2011B
    175860, 175863,
    175865, 175866, 175872, 175873, 175874, 175877, 175881, 175886, 175887,
    175888, 175906, 175910, 175921, 175973, 175974, 175975, 175976, 175990,
    176023, 176161, 176163, 176165, 176167, 176169, 176201, 176202, 176206,
    176207, 176286, 176289, 176304, 176308, 176309, 176467, 176468, 176469,
    176545, 176547, 176548, 176697, 176701, 176702, 176765, 176771, 176795,
    176796, 176797, 176799, 176801, 176805, 176807, 176841, 176844, 176848,
    176850, 176860, 176868, 176885, 176886, 176889, 176928, 176929, 176933,
    176959, 176982, 177053, 177074, 177088, 177095, 177096, 177131, 177138,
    177139, 177140, 177141, 177183, 177184, 177201, 177222, 177317, 177318,
    177319, 177449, 177452};
  */
  /*
    // **** 2.2/fb, Sep2nd *****
    // May10ReReco
    160431, 160577, 160578, 160871, 160872, 160873, 160874, 160939, 160940,
    160942, 160943, 160955, 160956, 160957, 160998, 161008, 161016, 161103,
    161106, 161107, 161113, 161116, 161117, 161119, 161156, 161176, 161217,
    161222, 161223, 161233, 161310, 161311, 161312, 162762, 162765, 162803,
    162808, 162811, 162822, 162825, 162826, 162828, 162909, 163046, 163069,
    163071, 163078, 163232, 163233, 163234, 163235, 163237, 163238, 163252,
    163255, 163261, 163270, 163286, 163289, 163296, 163297, 163300, 163301,
    163302, 163332, 163333, 163334, 163337, 163338, 163339, 163340, 163358,
    163369, 163370, 163371, 163372, 163374, 163375, 163376, 163378, 163385,
    163387, 163402, 163475, 163476, 163478, 163479, 163480, 163481, 163482,
    163483, 163582, 163583, 163584, 163585, 163586, 163587, 163588, 163589,
    163596, 163630, 163655, 163657, 163658, 163659, 163660, 163661, 163662,
    163663, 163664, 163668, 163738, 163757, 163758, 163759, 163760, 163761,
    163763, 163765, 163795, 163796, 163817, 163869,
    // PromptV4
    165088, 165098, 165099, 165102, 165103, 165120, 165121, 165205, 165208,
    165364, 165402, 165415, 165467, 165472, 165486, 165487, 165506, 165514,
    165548, 165558, 165567, 165570, 165617, 165620, 165633, 165970, 165993,
    166011, 166033, 166034, 166049, 166149, 166150, 166161, 166163, 166164,
    166346, 166374, 166380, 166408, 166429, 166438, 166462, 166486, 166502,
    166512, 166514, 166530, 166554, 166563, 166565, 166699, 166701, 166763,
    166781, 166782, 166784, 166787, 166839, 166841, 166842, 166859, 166860,
    166861, 166864, 166888, 166889, 166890, 166894, 166895, 166911, 166922,
    166923, 166946, 166950, 166960, 166966, 166967, 167039, 167041, 167043,
    167078, 167098, 167102, 167103, 167151, 167281, 167282, 167284, 167551,
    167673, 167674, 167675, 167676, 167740, 167746, 167754, 167784, 167786,
    167807, 167830, 167898, 167913,
    // Aug5ReReco
    170722, 170826, 170842, 170854, 170876, 170896, 170899, 170901, 171050,
    171091, 171098, 171102, 171106, 171117, 171156, 171178, 171219, 171274,
    171282, 171315, 171369, 171446, 171484, 171578, 171812, 171876, 171880,
    171895, 171897, 171921, 171926, 172014, 172024, 172033, 172163, 172208,
    172252, 172254, 172255, 172268, 172286, 172389, 172399, 172400, 172401,
    172411, 172478, 172619,
    // PromptV6
    172620, 172630, 172635, 172778, 172791, 172798, 172799, 172801, 172802,
    172819, 172822, 172824, 172847, 172865, 172868, 172949, 172951, 172952,
    172992, 172999, 173198, 173236, 173240, 173241, 173243, 173380, 173381,
    173389, 173406, 173430, 173438, 173439, 173657, 173658, 173659, 173660,
    173661, 173663, 173664, 173692};
  */
/*
// These are for the 1.1/fb file
    160431, 160577, 160578, 160871, 160872, 160873, 160874, 160939, 160940,
    160942, 160943, 160955, 160956, 160957, 160998, 161008, 161016, 161103,
    161106, 161107, 161113, 161116, 161117, 161119, 161156, 161176, 161217,
    161222, 161223, 161233, 161310, 161311, 161312, 162762, 162765, 162803,
    162808, 162811, 162822, 162825, 162826, 162828, 162909, 163046, 163069,
    163071, 163078, 163232, 163233, 163234, 163235, 163237, 163238, 163252,
    163255, 163261, 163270, 163286, 163289, 163296, 163297, 163300, 163301,
    163302, 163332, 163333, 163334, 163337, 163338, 163339, 163340, 163358,
    163369, 163370, 163371, 163372, 163374, 163375, 163376, 163378, 163385,
    163387, 163402, 163475, 163476, 163478, 163479, 163480, 163481, 163482,
    163483, 163582, 163583, 163584, 163585, 163586, 163587, 163588, 163589,
    163596, 163630, 163655, 163657, 163658, 163659, 163660, 163661, 163662,
    163663, 163664, 163668, 163738, 163757, 163758, 163759, 163760, 163761,
    163763, 163765, 163795, 163796, 163817, 163869, 165088, 165098, 165099,
    165102, 165103, 165120, 165121, 165205, 165208, 165364, 165402, 165415,
    165467, 165472, 165486, 165487, 165506, 165514, 165548, 165558, 165567,
    165570, 165617, 165620, 165633, 165970, 165993, 166011, 166033, 166034,
    166049, 166149, 166150, 166161, 166163, 166164, 166346, 166374, 166380,
    166408, 166429, 166438, 166462, 166486, 166502, 166512, 166514, 166530,
    166554, 166563, 166565, 166699, 166701, 166763, 166781, 166782, 166784,
    166787, 166839, 166841, 166842, 166859, 166860, 166861, 166864, 166888,
    166889, 166890, 166894, 166895, 166911, 166922, 166923, 166946, 166950,
    166960, 166966, 166967, 167039, 167041, 167043, 167078, 167098, 167102,
    167103, 167151, 167281, 167282, 167284, 167551, 167673, 167674, 167675,
    167740, 167746, 167754, 167784, 167786, 167807, 167830, 167898, 167913};
*/

  const int ngoodruns = sizeof(a_goodruns)/sizeof(a_goodruns[0]);
  set<int> goodruns;
  for (int i = 0; i != ngoodruns; ++i) {
    goodruns.insert(a_goodruns[i]);
  }
  set<pair<int, int> > nolums;

  if (true) {
    for (set<int>::const_iterator it = goodruns.begin(); it != goodruns.end();
	 ++it) {
      cout << *it << ", ";
    }
    cout << endl;
  }

  ifstream f(filename, ios::in);
  assert(f.is_open());
  float secLS = 2.3310e+01;
  string s;
  int rn, ls, ifoo;
  float del, rec, ffoo;
  assert(getline(f, s, '\r')); //assert(f >> s);
  cout << endl << "string: " << s << " !" << endl << flush;
  bool v1 = (s==string("run,ls,delivered,recorded"));
  bool v2 = (s==string("Run,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)"));
  assert(v1 || v2);
  int nls(0);
  double lumsum(0);
  double lumsum_good(0);
  double lumsum_json(0);
  bool skip(false);
  while ((v1 && f >> s &&
	  sscanf(s.c_str(),"%d,%d,%f,%f",&rn,&ls,&del,&rec)==4) ||
	 (v2 && getline(f, s, '\r') &&
	  (sscanf(s.c_str(),"%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,"
		  "%f,%f,%f", &rn,&ls,&ifoo, &ifoo,&ifoo,&ifoo,
		  &ifoo,&ifoo,&ifoo, &ffoo,&del,&rec)==12 ||
	   (skip=true)))) {

    // LS is not STABLE BEAMS but something else:
    // ADJUST, BEAM DUMP, FLAT TOP, INJECTION PHYSICS BEAM, N/A, RAMP DOWN,
    // SETUP, SQUEEZE
    if (skip) {
      skip = false;
      continue;
    }

    assert(_lums[rn][ls]==0);
    // Try to get this in units of pb-1
    // apparently it is given in s^-1 cm^-2
    //double lum = lvtx * secLS * 1e-36 ;
    //if (lum==0) lum = lhf * secLS * 1e-36 ;
    // lumiCalc.py returns lumi in units of mub-1 (=>nb-1=>pb-1) 
    double lum = rec*1e-6;
    //double lum2 = lhf * secLS * 1e-36 ;
    //if (lum2==0) lum2 = lvtx * secLS * 1e-36 ;
    double lum2 = del*1e-6;
    //assert(lum!=0);
    if (lum==0 && goodruns.find(rn)!=goodruns.end() &&
	_json[rn][ls]==1) {
      //cerr << "Warning: run " << rn << " LS " << ls << " is zero!" << endl;
      nolums.insert(make_pair<int, int>(rn,ls));
    }

    _lums[rn][ls] = lum;
    _lums2[rn][ls] = lum2;
    lumsum += lum;
    if (goodruns.find(rn)!=goodruns.end()) // Apr 17
      lumsum_good += lum;
    if (_json[rn][ls])
      lumsum_json += lum;
    ++nls;
    assert(nls<100000000);
  }

  cout << "Called loadLumi(\"" << filename << "\"):" << endl;
  cout << "Loaded " << _lums.size() << " runs with "
       << nls << " lumi sections containing "
       << lumsum << " pb-1 of data,\n of which "
       << lumsum_good << " pb-1 is in good runs ("
       << 100.*lumsum_good/lumsum << "%)"<< endl;
  cout << "This corresponds to " << nls*secLS/3600
       << " hours of data-taking" << endl;
  cout << "The JSON file contains "
       << lumsum_json << " pb-1 ("
       << 100.*lumsum_json/lumsum << "%)"<< endl;

  // Report any empty lumi section
  if (nolums.size()!=0) {
    cout << "Warning, found " << nolums.size() << " non-normalizable LS:";
    for (set<pair<int, int> >::const_iterator it = nolums.begin();
	 it != nolums.end(); ++it) {

      cout << " ["<<it->first<<","<<it->second;
      set<pair<int, int> >::const_iterator jt = it; ++jt;
      if (jt->first!=it->first || jt->second!=it->second+1)
	cout << "]";
      else {
	for (int i = 0; jt!=nolums.end() && jt->first==it->first
	       && jt->second==it->second+i+1; ++i, ++jt) {};
	it = --jt;
	cout << "-" << it->second << "]";
      }
    } // for it
    cout << endl;
  } // nolums
  
} // loadLumi


/*
void fillHistos::loadPrescale(const char *filename, string trigger) {

  ifstream f(filename, ios::in);
  assert(f.is_open());

  // Read in header information
  string s, trig, s1, s2, s3;
  f >> s >> trig; assert(s=="Trigger:");
  f >> s1 >> s2 >> s3;
  assert(s1=="Run");
  assert(s2=="L1_prescale");
  assert(s3=="HLT_prescale");

  // Read in data
  int rn, l1, hlt;
  while (f >> rn >> l1 >> hlt) {
    //assert(l1==1||l1==25);
    //assert(hlt>=1);
    _prescales[trigger][rn] = l1 * hlt;
  }

  cout << "Loaded prescales for " << trigger << " (" << trig << ") in "
       << _prescales[trigger].size() << " runs" << endl;

}
*/

void fillHistos::loadPUProfiles(const char *datafile, const char *mcfile) {

  cout << "Processing loadPUProfiles(\"" << datafile << "\")..." << endl;

  TDirectory *curdir = gDirectory;
  
  // Load pile-up files and hists from them
  TFile *fpudist = new TFile(datafile, "READ");
  assert(fpudist && !fpudist->IsZombie());
  TFile *fpumc = new TFile(mcfile,"READ");
  assert(fpumc && !fpumc->IsZombie());
  
  //pudist = (TH1D*)fpudist->Get("pileup"); assert(pudist);
  pumc = (TH1D*)fpumc->Get("hpu370"); assert(pumc);

  // Normalize
  //pudist->Scale(1./pudist->Integral());
  pumc->Scale(1./pumc->Integral());

  // For data, load each trigger separately
  //vector<string> triggers;
  _triggers.push_back("jt30");
  _triggers.push_back("jt60");
  //_triggers.push_back("jt80");
  _triggers.push_back("jt110");
  //_triggers.push_back("jt150");
  _triggers.push_back("jt190");
  _triggers.push_back("jt240");
  _triggers.push_back("jt300");
  _triggers.push_back("jt370");

  for (unsigned int itrg = 0 ; itrg != _triggers.size(); ++itrg) {

    const char *t = _triggers[itrg].c_str();
    pudist[t] = (TH1D*)fpudist->Get(Form("pileup_%s",t)); assert(pudist[t]);
    pudist[t]->Scale(1./pudist[t]->Integral());
  }

  curdir->cd();
} // loadPUProfiles

void fillHistos::loadECALveto(const char *file) {
  
  cout << "Processing loadECALveto(\"" << file << "\")..." << endl;

  TDirectory *curdir = gDirectory;
  
  TFile *fe = new TFile(file, "READ");
  assert(fe && !fe->IsZombie());
  
  ecalveto = (TH2F*)fe->Get("ecalveto"); assert(ecalveto);
  
  curdir->cd();
} // loadECALveto
