#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>
#include <assert.h>
#include <TGraph.h>
#include <THStack.h>
#include <TH1.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <ProfileProjector.h>
#include <tdrstyle_mod12.C>
#include <TAttMarker.h>
#include <sstream>
#include <HistScripts.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void runSingleHistos(string readMcFile="mc_histos.root" 
  , string filePath=
  "./")
{
  cout << "Processing runDrawHstos.C.." << endl;


  int puInspect = 0;
  int etaProfile = 0;
  int ptProfile = 0;
  int ptPrafile = 1;
  int ptDifference = 0;
  int debugger = 0;
  int angle = 0;
  // This indicates the file to be opened
  gROOT->Reset();
  gROOT->ProcessLine(".L QCDModules/QCDEvent.cc++");
  gROOT->ProcessLine(".L QCDModules/QCDJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDMET.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDPFJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDEventHdr.cc++");
  gROOT->ProcessLine(".L QCDModules/LorentzVector.h++");

  gROOT->ProcessLine(".L ProfileProjector.h+");

  string mcPath(filePath);
  mcPath += readMcFile;

  cout << mcPath << endl;
  ProfileProjector mcProcessor(mcPath,"mc");

  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("testcanvas","Generic Canvas",600,600);
  vector<TH1D*> mcPuHists = mcProcessor.GetPileUp();

  canvas->cd();
  setTDRStyle();
  canvas->UseCurrentStyle();
  if (puInspect){
    canvas->Divide(3,2,0,0);
    for (int i=0; i<6; i++){
      mcPuHists[i]->Scale(1./mcPuHists[i]->Integral());
      canvas->cd(i+1);
      mcPuHists[i]->Draw("");
      Double_t meanMC = mcPuHists[i]->GetMean();
      Double_t meanMCError = mcPuHists[i]->GetMeanError();
      std::stringstream textMC;
      textMC << "MC: " << meanMC << " #pm " << meanMCError;
      TLatex info;
      info.SetNDC();
      info.SetTextFont(43);
      info.SetTextSize(20);
      info.DrawLatex(0.5,0.6,textMC.str().c_str());
    }
  } else if (ptProfile){
    canvas->SetLogx();
    canvas->SetLogy();
    TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
    firstHisto->GetXaxis()->SetNoExponent();
    firstHisto->GetXaxis()->SetMoreLogLabels();
    firstHisto->Draw("");
    //pythiaFinal();
  } else if (ptPrafile){
    canvas->SetLogx();
    canvas->SetLogy();
    TH1D *firstHisto = mcProcessor.PtProfileScaled(1,0,-1);
    firstHisto->GetXaxis()->SetNoExponent();
    firstHisto->GetXaxis()->SetMoreLogLabels();
    firstHisto->Draw("");
    //pythiaFinal();
  } else if (etaProfile){
    TH1D *firstHisto = mcProcessor.PtProfile(0,0,-1);
    firstHisto->Draw("");
  } else if (ptDifference){
    canvas->SetLogx();
    canvas->SetLogy();
    TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
    firstHisto->Scale( 1./firstHisto->Integral()  );
    firstHisto->SetName("muu");
    firstHisto->Draw("same");
  } else if (debugger){
    THStack *tutkain1 = mcProcessor.EnergyFractions(1,0,-1,0,1);
    canvas->SetLogx();
    tutkain1->Draw("");
    mcProcessor.GetLeg()->Draw();
    canvas->Modified();
    canvas->SetSelected(canvas);
  } else if (angle){
    THStack *tutkain1 = mcProcessor.EnergyFractions(0,0,-1,0,1);
    canvas->cd();
    tutkain1->Draw("");
    mcProcessor.GetLeg()->Draw();
    canvas->Modified();
    canvas->SetSelected(canvas);
  }
}

