#ifndef HistScripts_h
#define HistScripts_h

#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>
#include <cassert>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

// Scale x-bins for logarithmic x-axis
void BinLogX(TH1 *hist) {
    TAxis *axis = hist->GetXaxis();
    int bins = axis->GetNbins();

    Axis_t from = TMath::Log10( axis->GetXmin() );
    Axis_t to = TMath::Log10( axis->GetXmax() );
    Axis_t width = (to - from) / bins;
    Axis_t *newBins = new Axis_t[bins + 1]; 

    for (int i = 0; i <= bins; i++){
        newBins[i] = TMath::Power(10, from + i * width);
    }   
    axis->Set(bins, newBins);
    delete newBins;
}

// A function that helps building a 1D-histogram
TH1 *HistBuilder( string title, Long64_t bins, Double_t lowLim, Double_t hiLim,
    string name, string xTitle, string yTitle, TCanvas *canvas, Int_t canvPos, 
    Int_t logY=1, Int_t logX=1, Double_t offSet = 1.2)
{
    TH1 *hist = new TH1D(name.c_str(),title.c_str(),bins,lowLim,hiLim);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetTitleOffset(offSet);
    if (logY){
        canvas->cd(canvPos)->SetLogy();
    }
    if (logX){
        canvas->cd(canvPos)->SetLogx();
        hist->GetXaxis()->SetMoreLogLabels();
        hist->GetXaxis()->SetNoExponent();
        BinLogX(hist);
    }
    return hist;
}

// A simpler version with no canvas-related actions
TH1 *HistBuilder( string title, Long64_t bins, Double_t lowLim, Double_t hiLim,
    string name, string xTitle, string yTitle, Int_t logX=1, Double_t offSet = 1.2)
{
    TH1 *hist = new TH1D(name.c_str(),title.c_str(),bins,lowLim,hiLim);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetTitleOffset(offSet);
    if (logX){
        hist->GetXaxis()->SetMoreLogLabels();
        hist->GetXaxis()->SetNoExponent();
        BinLogX(hist);
    }
    return hist;
}


#endif
