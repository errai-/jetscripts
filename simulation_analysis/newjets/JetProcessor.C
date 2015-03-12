#define JetProcessor_cxx
#include "JetProcessor.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "THStack.h"
#include "tdrstyle_mod14.C"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"

using std::cout;
using std::endl;


void stackModify(TH1D *setter){
  setter->GetXaxis()->SetTitle("p_{T} (GeV)");
  setter->GetYaxis()->SetTitle("Jet number fraction");
  setter->SetStats(0);
  setter->GetXaxis()->SetMoreLogLabels();
  setter->GetXaxis()->SetNoExponent();
}

void JetProcessor::Loop()
{
  if (fChain == 0) return;

  static const int ptBins = 48.;
  const double ptRange[ptBins+1]=
  //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
  //3637, 3832, 4037};//
  TH1D *gFrac = new TH1D("Gluon fractions","",ptBins,ptRange);
  TH1D *bFrac = new TH1D("b fractions","",ptBins,ptRange);
  TH1D *cFrac = new TH1D("c fractions","",ptBins,ptRange);
  TH1D *lFrac = new TH1D("light quark fractions","",ptBins,ptRange);
  TH1D *uFrac = new TH1D("unknown fractions","",ptBins,ptRange);

  size_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (size_t j=0; j!=nentries;j++) {
    Long64_t ientry = LoadTree(j);
    if (ientry < 0) break;
    nb = fChain->GetEntry(j);   nbytes += nb;
    
    size_t klim = fJets_;
    for (size_t k=0; k!=klim;k++){
      TLorentzVector P4;
      P4.SetPxPyPzE(fJets_fP4_fCoordinates_fX[k],fJets_fP4_fCoordinates_fY[k],
        fJets_fP4_fCoordinates_fZ[k],fJets_fP4_fCoordinates_fT[k]);
      if (fJets_fFlavour[k] == 21){
        gFrac->Fill( P4.Pt() );
      } else if (fJets_fFlavour[k] == 5){
        bFrac->Fill( P4.Pt() );
      } else if (fJets_fFlavour[k] == 4){
        cFrac->Fill( P4.Pt() );
      } else if (fJets_fFlavour[k] > 0){
        lFrac->Fill( P4.Pt() );
      } else {
        uFrac->Fill( P4.Pt() );
      }
    }
  }

  TH1D *sumFrac = new TH1D("All fractions","",ptBins,ptRange);
  sumFrac->Add(gFrac); sumFrac->Add(lFrac);  sumFrac->Add(bFrac); 
  sumFrac->Add(cFrac);  sumFrac->Add(uFrac);
  
  gFrac->Divide(sumFrac); bFrac->Divide(sumFrac); 
  cFrac->Divide(sumFrac);
  lFrac->Divide(sumFrac); uFrac->Divide(sumFrac);

  THStack* fracs = new THStack("flavour stacks","");

  TCanvas *canv = tdrCanvas("c1",gFrac,12,0,1);
  setTDRStyle();
  stackModify( gFrac );
  canv->UseCurrentStyle();

  gFrac->SetFillColor(kRed-4);
  bFrac->SetFillColor(kMagenta-4);
  cFrac->SetFillColor(kGreen-4);
  lFrac->SetFillColor(kTeal);
  uFrac->SetFillColor(kBlue-4);

  fracs->Add(gFrac); fracs->Add(lFrac); fracs->Add(cFrac); 
  fracs->Add(bFrac); fracs->Add(uFrac);
  fracs->SetHistogram( gFrac );
  fracs->Draw();
  TLegend * leg = tdrLeg(0.2,0.2,0.5,0.4);
  leg->AddEntry(uFrac, "unknown");
  leg->AddEntry(bFrac, "b");
  leg->AddEntry(cFrac, "c");
  leg->AddEntry(lFrac, "light quark");
  leg->AddEntry(gFrac, "gluons");
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(0.6*0.07);

  latex.DrawLatex(0.8,0.85,"Pythia 8");
  //CMS_lumi(canv,1,1);
  canv->Modified();
}
