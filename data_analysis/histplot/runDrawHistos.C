#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>

#include "TROOT.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "THStack.h"
#include "TH1.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TAttMarker.h"

#include "ProfileProjector.h"
#include "tdrstyle_mod12.C"

using std::string;
using std::vector;
using std::cout;
using std::endl;

void runDrawHistos(string readMcFile, string readDtFile)
{
    int puInspect = 0;
    int etaProfile = 0;
    int ptProfile = 0;
    int ptDifference = 0;
    int debugger = 0;
    int angle = 1;

    ProfileProjector mcProcessor(readMcFile,"mc");
    ProfileProjector dtProcessor(readDtFile,"dt"); 

    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("testcanvas","Generic Canvas",600,600);
    vector<TH1D*> mcPuHists = mcProcessor.GetPileUp();
    vector<TH1D*> dtPuHists = dtProcessor.GetPileUp();

    canvas->cd();
    setTDRStyle();
    canvas->UseCurrentStyle();
    if (puInspect) {
        canvas->Divide(3,2,0,0);
        for (int i=0; i<6; i++) {
            mcPuHists[i]->Scale(1./mcPuHists[i]->Integral());
            dtPuHists[i]->Sumw2();
            dtPuHists[i]->Scale(1./dtPuHists[i]->Integral());
            canvas->cd(i+1);
            std::stringstream tmpStr("");
            tmpStr << "trigger" << i+1;
            mcPuHists[i]->SetTitle(tmpStr.str().c_str());
            mcPuHists[i]->Draw("");
            dtPuHists[i]->SetMarkerStyle(kCircle);
            dtPuHists[i]->SetMarkerColor(kRed);
            dtPuHists[i]->Draw("SAMEp");
            cout << mcPuHists[i]->GetTitle() << " " << dtPuHists[i]->GetTitle() << endl;
            Double_t meanMC = mcPuHists[i]->GetMean();
            Double_t meanMCError = mcPuHists[i]->GetMeanError();
            Double_t meanDT = dtPuHists[i]->GetMean();
            Double_t meanDTError = dtPuHists[i]->GetMeanError();
            std::stringstream textMC;
            std::stringstream textDT;
            textMC << "MC: " << meanMC << " #pm " << meanMCError;
            textDT << "DT: " << meanDT << " #pm " << meanDTError;
            TLatex info;
            info.SetNDC();
            info.SetTextFont(43);
            info.SetTextSize(20);
            info.DrawLatex(0.5,0.6,textMC.str().c_str());
            info.DrawLatex(0.5,0.5,textDT.str().c_str());
        }
        canvas->cd(0);
        //cmsFinal();
    } else if (ptProfile) {
        canvas->SetLogx();
        canvas->SetLogy();
        TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
        firstHisto->Draw("");
        TH1D* tmpHisto = dtProcessor.PtProfile(1,0,-1);
        tmpHisto->SetMarkerStyle(kCircle);
        tmpHisto->SetMarkerColor(kRed);
        tmpHisto->Draw("samep");
    } else if (etaProfile) {
        TH1D *firstHisto = mcProcessor.PtProfile(0,0,-1);
        firstHisto->Draw("");
        TH1D* tmpHisto = mcProcessor.PtProfile(0,1000,2000);
        tmpHisto->SetMarkerStyle(kCircle);
        tmpHisto->SetMarkerColor(kRed);
        tmpHisto->Draw("samep");
    } else if (ptDifference) {
        canvas->SetLogx();
        canvas->SetLogy();
        TH1D* tmpHisto = dtProcessor.PtProfile(1,0,-1);
        tmpHisto->Scale(1./tmpHisto->Integral());
        tmpHisto->SetMarkerStyle(kCircle);
        tmpHisto->SetMarkerColor(kRed);
        tmpHisto->Draw("p");
        TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
        firstHisto->Scale( 1./firstHisto->Integral()  );
        firstHisto->SetName("muu");
        firstHisto->Draw("same");
    } else if (debugger) {
        THStack *tutkain1 = mcProcessor.EnergyFractions(1,0,-1,0,1);
        THStack *tutkain2 = dtProcessor.EnergyFractions(1,0,-1,0,1);
        canvas->SetLogx();
        tutkain2->Draw("");
        tutkain1->Draw("samep");
        dtProcessor.GetLeg()->Draw();
        cmsFinal();
        canvas->Modified();
        canvas->SetSelected(canvas);
    } else if (angle) {
        THStack *tutkain1 = mcProcessor.EnergyFractions(0,0,-1,0,1);
        THStack *tutkain2 = dtProcessor.EnergyFractions(0,0,-1,1,1);
        canvas->cd();
        tutkain1->Draw("");
        tutkain2->Draw("samee1");
        mcProcessor.GetLeg()->Draw();
        canvas->Modified();
        canvas->SetSelected(canvas);
    }
}

