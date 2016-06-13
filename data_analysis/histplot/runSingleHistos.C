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

void runSingleHistos(string readFile="mc_histos.root")
{
    int puInspect = 0;
    int etaProfile = 0;
    int ptProfile = 0;
    int ptPrafile = 0;
    int ptDifference = 0;
    int debugger = 1;
    int angle = 1;

    ProfileProjector mcProcessor(readFile,"mc");

    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("testcanvas","Generic Canvas",600,600);
    vector<TH1D*> mcPuHists = mcProcessor.GetPileUp();

    canvas->cd();
    setTDRStyle();
    canvas->UseCurrentStyle();
    if (puInspect) {
        canvas->Divide(3,2,0,0);
        for (int i=0; i<6; i++) {
            mcPuHists[i]->Scale(1./mcPuHists[i]->Integral());
            canvas->cd(i+1);
            mcPuHists[i]->Draw("");
            Double_t meanMC = mcPuHists[i]->GetMean();
            Double_t meanMCError = mcPuHists[i]->GetMeanError();
            std::stringstream textMC;
            textMC << "DT: " << meanMC << " #pm " << meanMCError;
            TLatex info;
            info.SetNDC();
            info.SetTextFont(43);
            info.SetTextSize(20);
            info.DrawLatex(0.5,0.6,textMC.str().c_str());
        }
    } else if (ptProfile) {
        canvas->SetLogx();
        canvas->SetLogy();
        TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
        firstHisto->GetXaxis()->SetNoExponent();
        firstHisto->GetXaxis()->SetMoreLogLabels();
        firstHisto->Draw("");
        //pythiaFinal();
    } else if (ptPrafile) {
        canvas->SetLogx();
        canvas->SetLogy();
        TH1D *firstHisto = mcProcessor.PtProfileScaled(1,0,-1);
        firstHisto->GetXaxis()->SetNoExponent();
        firstHisto->GetXaxis()->SetMoreLogLabels();
        firstHisto->Draw("HIST");
        //pythiaFinal();
    } else if (etaProfile) {
        TH1D *firstHisto = mcProcessor.PtProfile(0,0,-1);
        firstHisto->Draw("");
    } else if (ptDifference) {
        canvas->SetLogx();
        canvas->SetLogy();
        TH1D *firstHisto = mcProcessor.PtProfile(1,0,-1);
        firstHisto->Scale( 1./firstHisto->Integral()  );
        firstHisto->SetName("muu");
        firstHisto->Draw("same");
    } else if (debugger) {
        THStack *stackPlot = mcProcessor.EnergyFractions(1,0,-1,0,1);
        canvas->SetLogx();
        stackPlot->Draw("");
        mcProcessor.GetLeg()->Draw();
        canvas->Modified();
        canvas->SetSelected(canvas);
    } else if (angle) {
        THStack *stackPlot = mcProcessor.EnergyFractions(0,0,-1,0,1);
        canvas->cd();
        stackPlot->Draw("");
        mcProcessor.GetLeg()->Draw();
        canvas->Modified();
        canvas->SetSelected(canvas);
    }
}

