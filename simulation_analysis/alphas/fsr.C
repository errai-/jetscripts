// Check FSR difference between Pythia6 and Pythia8
// input files from pure generator level MC from Hannu
#include "TFile.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "tdrstyle_mod14.C"

#include <string>
#include <map>

using namespace std;

void fsr(int ins = -1, double pt = 200., std::string sample = "dijet") {

  setTDRStyle();

  TDirectory *curdir = gDirectory;

  const char *cs = sample.c_str();
  map<string, const char *> title;
  title["dijet"] = "Dijet";
  title["gamjet"] = "#gamma+jet";

  TFile *fp8 = new TFile("alphafracs_p8.root","READ");
  assert(fp8 && !fp8->IsZombie());

  TFile *fp6 = new TFile("alphafracs_p6.root","READ");
  assert(fp6 && !fp6->IsZombie());

  TFile *fhw = new TFile("alphafracs_hwpp.root","READ");
  assert(fhw && !fhw->IsZombie());

  curdir->cd();

  TProfile *p8_10 = (TProfile*)fp8->Get("prof10"); assert(p8_10);
  TProfile *p8_15 = (TProfile*)fp8->Get("prof15"); assert(p8_15);
  TProfile *p8_20 = (TProfile*)fp8->Get("prof20"); assert(p8_20);
  TProfile *p8_30 = (TProfile*)fp8->Get("prof30"); assert(p8_30);
  TH1D *h8_0 = p8_10->ProjectionX("h8_10");

  TProfile *p6_10 = (TProfile*)fp6->Get("prof10"); assert(p6_10);
  TProfile *p6_15 = (TProfile*)fp6->Get("prof15"); assert(p6_15);
  TProfile *p6_20 = (TProfile*)fp6->Get("prof20"); assert(p6_20);
  TProfile *p6_30 = (TProfile*)fp6->Get("prof30"); assert(p6_30);
  TH1D *h6_0 = p6_10->ProjectionX("h6_10");

  TProfile *hw_10 = (TProfile*)fhw->Get("prof10"); assert(hw_10);
  TProfile *hw_15 = (TProfile*)fhw->Get("prof15"); assert(hw_15);
  TProfile *hw_20 = (TProfile*)fhw->Get("prof20"); assert(hw_20);
  TProfile *hw_30 = (TProfile*)fhw->Get("prof30"); assert(hw_30);
  TH1D *hw_0 = hw_10->ProjectionX("hw_10");

  TH1D *dt_0 = (TH1D*)hw_0->Clone("dt_0");
  TH1D *dt_10 = hw_10->ProjectionX("dt_10");
  TH1D *dt_15 = hw_15->ProjectionX("dt_15");
  TH1D *dt_20 = hw_20->ProjectionX("dt_20");
  TH1D *dt_30 = hw_30->ProjectionX("dt_30");

  const int ns = 3;//4;
  const int np = 5;
  TH1D* ps[ns][np] = {{h8_0, p8_10, p8_15, p8_20, p8_30},
		      {h6_0, p6_10, p6_15, p6_20, p6_30},
		      {hw_0, hw_10, hw_15, hw_20, hw_30}};//,
  //{dt_0, dt_10, dt_15, dt_20, dt_30}};
  double alpha[np] = {0, 0.10, 0.15, 0.20, 0.30};
  int markers[ns][2] = {{kFullSquare, kFullCircle},
			{kOpenSquare, kOpenCircle},
			{kOpenDiamond, kOpenStar}};//,
  //{kDot, kDot}};
  int colors[np] = {kBlack, kRed, kOrange+2, kGreen+2, kBlue};

  assert(ins>=-1 && ins<ns);

  // Approximate data as 1:1 mixture of P6 and Herwig++
  /*
  double whw = 0.5;
  for (int j = 0; j != np; ++j) {
    for (int k = 1; k != hw_0->GetNbinsX()+1; ++k) {
      ps[ns-1][j]->SetBinContent(k, (1-whw)*ps[1][j]->GetBinContent(k) +
				 whw*ps[2][j]->GetBinContent(k));
      ps[ns-1][j]->SetBinError(k, (1-whw)*ps[1][j]->GetBinError(k) +
			       whw*ps[2][j]->GetBinError(k));
    }
  }
  */

  TGraphErrors *gas[ns];
  TGraphErrors *ga = new TGraphErrors(4);
  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,0.35);
  
  for (int i = 0; i != ns; ++i) {
    for (int k = 1; k != p8_10->GetNbinsX()+1; ++k) {
      for (int j = 1; j != np; ++j) {
	ga->SetPoint(j, alpha[j], ps[i][j]->GetBinContent(k));
	ga->SetPointError(j, 0., ps[i][j]->GetBinError(k));
      }
      ga->Fit(f1, "QRN");
      ps[i][0]->SetBinContent(k, f1->GetParameter(0));
      ps[i][0]->SetBinError(k, f1->GetParError(0));

      if (p8_10->FindBin(pt)==k) {
	gas[i] = (TGraphErrors*)ga->Clone(Form("ga_%d",i));
      }
    }
  }

  TH1D *h = new TH1D("h",";p_{T,parton} (GeV);"
		     "#LTp_{T,gen} / p_{T,parton}#GT", 100, 0, 900);
  h->SetMaximum(1.06);//1.03);
  h->SetMinimum(0.91);//0.96);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  extraText = "Simulation";
  //extraText2 = "Preliminary";
  if (ins==0) lumi_13TeV = "Pythia8";
  if (ins==1) lumi_13TeV = "Pythia6";
  if (ins==2) lumi_13TeV = "Herwig++";
  if (ins==-1) lumi_13TeV = "Herwig++ / Pythia8 / Pythia6";

  TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);
  tex->DrawLatex(0.18,0.87,title[cs]);
  tex->DrawLatex(0.18,0.80,"Anti-k_{T} R=0.5");
  tex->DrawLatex(0.18,0.75,"|#eta| < 1.3");
  tex->DrawLatex(0.0,0.01,"#copyright Hannu Siikonen");

  for (int i = 0; i != ns; ++i) {
    for (int j = 0; j != np; ++j) {
      if ((ins == -1 && (j==0 || j==np-1)) || i == ins)
	tdrDraw(ps[i][j], "P", markers[i][j==0 ? 0 : 1], colors[j]);
    }
  }
  //tdrDraw(h8_0,"P",kFullSquare,kBlack);
  //tdrDraw(p8_10,"P",kFullCircle,kRed);
  //tdrDraw(p8_15,"P",kFullCircle,kOrange+2);
  //tdrDraw(p8_20,"P",kFullCircle,kGreen+1);
  //tdrDraw(p8_30,"P",kFullCircle,kBlue);

  TF1 *f2 = new TF1("f2","[0] + [1]*pow(x,[2])",100,840);
  f2->SetParameters(1,-0.1,-0.5);
  for (int i = 0; i != ns; ++i) {
    for (int j = 0; j != np; ++j) {
      ps[i][j]->Fit(f2,"QRN");
      f2->SetLineColor(ps[i][j]->GetMarkerColor());
      if ((ins == -1 && (j==0 || j==np-1)) || i == ins)
	f2->DrawClone("SAME");
    }
  }

  if (ins!=-1) {
    TLegend *leg = tdrLeg(0.70,0.65,0.90,0.90);
    leg->AddEntry(ps[ins][0],"#alpha_{max}#rightarrow0","PL");
    leg->AddEntry(ps[ins][1],"#alpha<0.10","PL");
    leg->AddEntry(ps[ins][2],"#alpha<0.15","PL");
    leg->AddEntry(ps[ins][3],"#alpha<0.20","PL");
    leg->AddEntry(ps[ins][4],"#alpha<0.30","PL");
  }
  if (ins==-1) {
    TLegend *ll = tdrLeg(0.50,0.70,0.70,0.90);
    ll->SetHeader("#alpha#rightarrow0");
    ll->AddEntry(ps[0][0],"","PL");
    ll->AddEntry(ps[1][0],"","PL");
    ll->AddEntry(ps[2][0],"","PL");
    //ll->AddEntry(ps[3][0],"","PL");
    TLegend *lr = tdrLeg(0.60,0.70,0.80,0.90);
    lr->SetHeader("#alpha<0.30");
    lr->AddEntry(ps[0][np-1],"   Pythia 8","PL");
    lr->AddEntry(ps[1][np-1],"   Pythia 6","PL");
    lr->AddEntry(ps[2][np-1],"   Herwig++","PL");
    //lr->AddEntry(ps[3][np-1],"   \"DATA\"","PL");
  }

  c1->SaveAs(Form("pdf/fsr_ins%d_%s.pdf",ins,cs));


  TH1D *h2 = new TH1D("h2",";#alpha_{max};#LTp_{T,gen} / p_{T,parton}#GT",
		      10,0,0.4);
  h2->SetMinimum(0.91);
  h2->SetMaximum(1.05);

  lumi_13TeV = "Herwig++ / Pythia8 / Pythia6";  
  TCanvas *c2 = tdrCanvas("c2",h2,2,0,kSquare);
  tex->DrawLatex(0.18,0.87,title[cs]);
  tex->DrawLatex(0.18,0.80,"Anti-k_{T} R=0.5");
  tex->DrawLatex(0.18,0.75,"|#eta| < 1.3");
  tex->DrawLatex(0.0,0.01,"#copyright Hannu Siikonen");

  TF1 *f3 = new TF1("f3","[0]+[1]*x+[2]*x*x",0,0.35);
  f3->SetLineStyle(kDashed);

  for (int i = 0; i != ns; ++i) {
    
    tdrDraw(gas[i], "P", markers[i][1], colors[i]);

    gas[i]->Fit(f1,"QRN");
    f1->SetLineColor(gas[i]->GetLineColor());
    f1->DrawClone("SAME");

    gas[i]->Fit(f3,"QRN");
    f3->SetLineColor(gas[i]->GetLineColor());
    f3->DrawClone("SAME");
  }

  int ipt = dt_0->FindBin(pt);
  double ptmin = dt_0->GetBinLowEdge(ipt);
  double ptmax = dt_0->GetBinLowEdge(ipt+1);
  tex->DrawLatex(0.18, 0.18, Form("%1.0f<p_{T}<%1.0f GeV",
				  ptmin, ptmax));

  TLegend *leg = tdrLeg(0.70,0.70,0.90,0.90);
  leg->AddEntry(gas[0], "Pythia 8", "PL");
  leg->AddEntry(gas[1], "Pythia 6", "PL");
  leg->AddEntry(gas[2], "Herwig++", "PL");
  //leg->AddEntry(gas[3], "\"DATA\"", "PL");

  c2->SaveAs(Form("pdf/fsr_vsalpha_%s.pdf",cs));
}
