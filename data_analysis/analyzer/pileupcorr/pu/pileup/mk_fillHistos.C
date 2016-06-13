// Purpose: Fill inclusive jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: August 15, 2011
{

  gROOT->ProcessLine(".exception");

  // compile code
  gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L runHistos.C+");
  gROOT->ProcessLine(".L basicHistos.C+");

  // For JEC residual (and pile-up)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");
  // For new JEC uncertainty
  //gROOT->ProcessLine(".L jecsys/ErrorTypes.cpp+");
  //gROOT->ProcessLine(".L jecsys/L3Corr.cpp+");

  gROOT->ProcessLine(".L fillHistos.C++g");

  gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  cafe::Config *cfg = new cafe::Config("pfjet");
  string type = cfg->get("type","DATA");
  delete cfg;

  #include "settings.h"
  #include <string>

  std::string algo = (_algo=="AK7" ? "ak7" : "ak5");
  const char *a = algo.c_str();

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",a));
  if (type=="DATA") {
    // Retrieved with rfcp /castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/July8th/InclusiveJetsTree_data.root InclusiveJetsTree_July8th.root
    //c->AddFile("data/InclusiveJetsTree_July8th.root/ak5/ProcessedTree");//1/fb
    //c->AddFile("data2/InclusiveJetsTree_Sep2nd.root/ak5/ProcessedTree");//2/fb
    //c->AddFile("data2/InclusiveJetsTree_Oct11th_Run2011A_May10.root/ak5/ProcessedTree");//3/fb
    //c->AddFile("data2/InclusiveJetsTree_Oct11th_Run2011A_Aug05.root/ak5/ProcessedTree");//3/fb
    //c->AddFile("data2/InclusiveJetsTree_Oct11th_Run2011A_PromptV4.root/ak5/ProcessedTree");//3/fb
    //c->AddFile("data2/InclusiveJetsTree_Oct11th_Run2011A_PromptV6.root/ak5/ProcessedTree");//3/fb
    //c->AddFile("data2/InclusiveJetsTree_Oct11th_Run2011B_PromptV1.root/ak5/ProcessedTree");//3/fb
    //c->AddFile("data2/InclusiveJetsTree_Nov13th_Run2011A_May10.root/ak5/ProcessedTree");//5/fb
    //c->AddFile("data2/InclusiveJetsTree_Nov13th_Run2011A_Aug05.root/ak5/ProcessedTree");//5/fb
    //c->AddFile("data2/InclusiveJetsTree_Nov13th_Run2011A_PromptV4.root/ak5/ProcessedTree");//5/fb
    //c->AddFile("data2/InclusiveJetsTree_Nov13th_Run2011A_PromptV6.root/ak5/ProcessedTree");//5/fb
    //c->AddFile("data2/InclusiveJetsTree_Nov13th_Run2011B_PromptV1.root/ak5/ProcessedTree");//5/fb
    /*
    c->AddFile(Form("data2/InclusiveJetsTree_Jan16th_Run2011A_May10.root/%s/ProcessedTree",a));//5/fb, new JEC, JSON
    c->AddFile(Form("data2/InclusiveJetsTree_Jan16th_Run2011A_Aug05.root/%s/ProcessedTree",a));//5/fb, new JEC, JSON
    c->AddFile(Form("data2/InclusiveJetsTree_Jan16th_Run2011A_PromptV4.root/%s/ProcessedTree",a));//5/fb, new JEC, JSON
    c->AddFile(Form("data2/InclusiveJetsTree_Jan16th_Run2011A_PromptV6.root/%s/ProcessedTree",a));//5/fb, new JEC, JSON
    c->AddFile(Form("data2/InclusiveJetsTree_Jan16th_Run2011B_PromptV1.root/%s/ProcessedTree",a));//5/fb, new JEC, JSON
    */
    c->AddFile("data2/InclusiveJetsTree_Jan16th_Run2011A_May10.root");
    c->AddFile("data2/InclusiveJetsTree_Jan16th_Run2011A_Aug05.root");
    c->AddFile("data2/InclusiveJetsTree_Jan16th_Run2011A_PromptV4.root");
    c->AddFile("data2/InclusiveJetsTree_Jan16th_Run2011A_PromptV6.root");
    c->AddFile("data2/InclusiveJetsTree_Jan16th_Run2011B_PromptV1.root");
  }
  if (type=="MC") {
    if (_pthatbins) {
      cout << "Running over pthat bins" << endl;
      c->AddFile("data2/Pythia6/Pythia6-5to15_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-15to30_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-30to50_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-50to80_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-80to120_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-120to170_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-170to300_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-300to470_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-470to600_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-600to800_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-800to1000_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-1000to1400_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-1400to1800_ProcessedTree_mc.root");
      c->AddFile("data2/Pythia6/Pythia6-1800to3500_ProcessedTree_mc.root");
    }
    else {
      cout << "Running over flat sample" << endl;
      //c->AddFile("data/Summer11PythiaZ2Flat_InclusiveJetsTree_mc.root/ak5/ProcessedTree");
      //c->AddFile(Form("data2/Pythia6Flat_ProcessedTree_mc.root/%s/ProcessedTree",a));
      c->AddFile("data2/Pythia6Flat_ProcessedTree_mc.root");
    }
  }
  if (type=="HW") {
    //c->AddFile("data/Summer11Herwig23Flat_InclusiveJetsTree_mc.root/ak5/ProcessedTree");
    //c->AddFile(Form("data2/HerwigFlat_ProcessedTree_mc.root/%s/ProcessedTree",a));
    c->AddFile("data2/HerwigFlat_ProcessedTree_mc.root");
  }
  fillHistos *h = new fillHistos(c, type);
  h->Loop();
  delete h;

  // Trying to process all one after the other seems to crash and mess output
  // use multifill.sh instead
  /*
  // Process first MC
  TChain *cm = new TChain("ak5/ProcessedTree");
  cm->AddFile("data/Summer11PythiaZ2Flat_InclusiveJetsTree_mc.root/ak5/ProcessedTree");
  fillHistos *hm = new fillHistos(cm, "MC");
  hm->Loop();
  delete hm;

  // ...then data
  TChain *cd = new TChain("ak5/ProcessedTree");
  // Retrieved with rfcp /castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/July8th/InclusiveJetsTree_data.root InclusiveJetsTree_July8th.root
  cd->AddFile("data/InclusiveJetsTree_July8th.root/ak5/ProcessedTree"); //1/fb
  fillHistos *hd = new fillHistos(cd, "DATA");
  hd->Loop();
  delete hd;

  // ...then another MC
  TChain *ch = new TChain("ak5/ProcessedTree");
  ch->AddFile("data/Summer11Herwig23Flat_InclusiveJetsTree_mc.root/ak5/ProcessedTree");
  fillHistos *hh = new fillHistos(ch, "HW");
  hh->Loop();
  delete hh;
  */
}
