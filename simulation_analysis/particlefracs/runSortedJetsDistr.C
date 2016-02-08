// This file uses the sortedjets data and studies the quark and gluon jets

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// ROOT, for histogramming.
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TPaveStats.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod14_pythia8.C"

int runSortedJetsDistr(std::string fileName)
{
    int choose = 8;

    // Create file on which histogram(s) can be saved.
    TFile *inFile = new TFile(fileName.c_str(), "READ");

    // Read histograms.
    vector<TH1D*> fractionHists;
    fractionHists.push_back( (TH1D*) inFile->Get( "gluonjet amount" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "quarkjet amount" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "gluonjet charge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "quarkjet charge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "gluonjet wcharge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "quarkjet wcharge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "gluonjet w2charge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "quarkjet w2charge" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "gluonjet w" ) );
    fractionHists.push_back( (TH1D*) inFile->Get( "quarkjet w" ) );

    TCanvas *canv = new TCanvas("c1","c1",600,600);
    setTDRStyle();
    canv->SetSelected(canv);
    canv->cd();
    // Show histogram
    fractionHists[choose]->Scale(1./fractionHists[choose]->Integral());
    fractionHists[choose+1]->Scale(1./fractionHists[choose+1]->Integral());
    if (choose == 0){
        fractionHists[0]->GetXaxis()->SetTitle("Particles in a jet");
        fractionHists[0]->GetYaxis()->SetTitle("Scaled distr.");
        fractionHists[0]->SetMaximum(0.055);
        fractionHists[1]->SetMaximum(0.055);
    } else if (choose == 2){
        fractionHists[2]->GetXaxis()->SetTitle("Absolute jet charge (e)");
        fractionHists[2]->GetYaxis()->SetTitle("Scaled distr.");
        fractionHists[2]->SetMaximum(0.25);
        fractionHists[3]->SetMaximum(0.25);
    } else if (choose == 4){
        fractionHists[4]->GetXaxis()->SetTitle("Weighted jet charge (pTi/pT)");
        fractionHists[4]->GetYaxis()->SetTitle("Scaled distr.");
        fractionHists[4]->SetMaximum(0.04);
        fractionHists[5]->SetMaximum(0.04);
    } else if (choose == 6){
        fractionHists[6]->GetXaxis()->SetTitle("Weighted jet charge (pTi/pT)^2");
        fractionHists[6]->GetYaxis()->SetTitle("Scaled distr.");
        fractionHists[6]->SetMaximum(0.15);
        fractionHists[7]->SetMaximum(0.15);
    } else if (choose == 8){
        fractionHists[8]->GetXaxis()->SetTitle("Squared weights (pTi/pT)^2");
        fractionHists[8]->GetYaxis()->SetTitle("Scaled distr.");
        fractionHists[8]->SetMaximum(0.1);
        fractionHists[9]->SetMaximum(0.1);
    }
    fractionHists[choose]->GetYaxis()->SetTitleOffset(1.85);
    fractionHists[choose]->GetXaxis()->SetTitleOffset(1.1);
    canv->UseCurrentStyle();
    fractionHists[choose+1]->SetMarkerStyle(kCircle);
    fractionHists[choose+1]->SetMarkerColor(kRed);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad1->Draw();
    pad1->cd();
    fractionHists[choose]->Draw();
    pad1->Update();
    TPaveStats *ps1 = (TPaveStats*)fractionHists[choose]->GetListOfFunctions()->FindObject("stats");
    if (choose == 0){
        ps1->SetX1NDC(0.6);
        ps1->SetX2NDC(0.8);
        ps1->SetY1NDC(0.5);
        ps1->SetY2NDC(0.7);
    } else {
        ps1->SetX1NDC(0.2);
        ps1->SetX2NDC(0.4);
        ps1->SetY1NDC(0.7);
        ps1->SetY2NDC(0.9);
    }
    ps1->SetTextColor(kBlue);
    pad1->Modified();
    canv->cd();
    pad2->Draw();
    pad2->cd();
    fractionHists[choose+1]->Draw("SAMEp");
    pad2->Update();
    TPaveStats *ps2 = (TPaveStats*)fractionHists[choose+1]->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6);
    ps2->SetX2NDC(0.8);
    ps2->SetY1NDC(0.7);
    ps2->SetY2NDC(0.9);
    ps2->SetTextColor(kRed);
    canv->cd();
    pad2->Modified();
    //pythiaFinal();
    std::cout << "\nDouble click on the histogram window to quit.\n";
    gPad->WaitPrimitive();

    // Done.
    return 0;
}
