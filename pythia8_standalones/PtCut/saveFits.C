#define pfhadrons_cxx
//#include "pfhadrons.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TLegend.h"
#include "TFile.h"

//#include "tdrstyle_mod12.C";
//#include "tdrstyle_mod14.C"
#include "TStyle.h"
#include "TGraphErrors.h"

using namespace std;

// energy below which ecal deposit is considered MIP, i.e.
// transition from eH to EH hadrons
const double mip_ecal = 0.8;//1.0;
// thresholds in ECAL and HCAL
const double minpt_hcal = 0.8;
const double minpt_ecal = 0.3;
const double minpt_corr = 0.8;


// H+E = f(T; H, E);
// x - variables, p - parameters
Double_t func_eh(Double_t *x, Double_t *p) {

  double T = x[0];
  double E = p[0];
  double H = p[1];
  // V2 with r+/-dr
  Double_t *pa = &p[6];
  Double_t *pb = &p[2];
  double a = max(0.33, (pa[0]+pa[1]*pow(T+pa[3],pa[2])) * max(0., T+pa[3])/T);
  double b = max(0.33, (pb[0]+pb[1]*pow(T+pb[3],pb[2])) * max(0., T+pb[3])/T);

  double EH = a*(T - H/b) + b*(T - E/a);

  return EH;
} // func_he

Double_t func(Double_t *x, Double_t *p) {
  
  double ht = x[0];
  double et = x[1];
  double n = p[0];
  double rh = p[1];
  double sh = p[2];
  double re = p[3];
  double se = sh; // assume same for now
  double rf = p[4];
  double sf = p[5];

  // Find nearest point to the ht = rh*fh, et = re*(1-fh) line
  // ie. et = re*(1-ht/rh) or ht = rh*(1-et/re)
  double ht0 = ht;
  double et0 = re*(1-ht0/rh);
  double et1 = et;
  double ht1 = rh*(1-et1/re);
  // Length of line segment 01
  double l = sqrt(pow(ht1-ht0,2) + pow(et1-et0,2));
  // Project vector 0x to line segment 01 to get length of projection
  double h = (ht-ht0)*(ht1-ht0) + (et-et0)*(et1-et0);
  // Move distance h from P0 along line segment l to get to nearest point
  double ht2 = h/l*(ht1-ht0)+ht0;
  double et2 = h/l*(et1-et0)+et0;

  double dh = ht - ht2;
  double de = et - et2;
  double df = (et2/re)/((ht2/rh)+(et2/re)) - rf;
   
  return (n*TMath::Gaus(dh,0,sh,kTRUE)*TMath::Gaus(de,0,se,kTRUE)
	  *TMath::Gaus(df,0,sf,kTRUE));
}

void saveFits(double pt = -5) {

  TDirectory *curdir = gDirectory;
  TFile *f = new TFile("pfhadrons.root","READ");
  assert(f && !f->IsZombie());

  TH3D *h3 = (TH3D*)f->Get("e3r"); assert(h3);

  // EH-hadrons response
  TGraphErrors *grh = new TGraphErrors(0);
  TGraphErrors *gre = new TGraphErrors(0);
  TGraphErrors *grs = new TGraphErrors(0);
  TGraphErrors *grr = new TGraphErrors(0);
  TGraphErrors *grd = new TGraphErrors(0);
  // EH-hadron ECAL fractions
  TGraphErrors *gfm = new TGraphErrors(0);
  TGraphErrors *gfs = new TGraphErrors(0);

  // H-hadrons
  TGraphErrors *ghr = new TGraphErrors(0);
  TGraphErrors *ghs = new TGraphErrors(0);
  TGraphErrors *gh4 = new TGraphErrors(0); // Shifted up by x5

  // E-hadrons
  TGraphErrors *ger = new TGraphErrors(0);
  TGraphErrors *ges = new TGraphErrors(0);

  // All hadrons
  TGraphErrors *gar = new TGraphErrors(0);
  TGraphErrors *gas = new TGraphErrors(0);


  TF2 *f2(0), *f2b(0);
  TF1 *f1(0); 
  double norm_mip = 0.045; // 0.065;
  double ptmin = (pt>0 ? min(pt,20.) : 3.);//8.);//20.);
  int ipt0 = (pt>0 ? h3->GetXaxis()->FindBin(pt) : 0);
  int iptmin = (pt>0 ? ipt0 : h3->GetXaxis()->FindBin(ptmin));
  double ipt1 = h3->GetXaxis()->FindBin(199.);
  bool goodfit = false;
  for (int ipt = h3->GetNbinsZ(); ipt != iptmin-1; --ipt) {

    if (ipt0!=0 && ipt!=ipt0) continue;
    double pt1 = h3->GetXaxis()->GetBinLowEdge(ipt);
    double pt2 = h3->GetXaxis()->GetBinLowEdge(ipt+1);
    double pt = 0.5*(pt1+pt2);
    double dpt = 0.5*(pt2-pt1);

    if (pt2>200) continue; // VLE+ME sample

    TH2D *h2 = (TH2D*)h3->Project3D("yz");
    h2->GetYaxis()->SetRangeUser(0,1.5);

    TProfile *px = h2->ProfileX();
    TProfile *py = h2->ProfileY();
    
    double minH0 = minpt_hcal/pt1;
    int iminH0 = px->FindBin(minH0);
    double minH = px->GetBinLowEdge(iminH0+1);
    double minE0 = mip_ecal/pt1;
    int iminE0 = py->FindBin(minE0);
    double minE = py->GetBinLowEdge(iminE0+1);

    // Select just one pT bin
    for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	      // normalize E=0 bin, which is due to MIPS
	      double ET = h2->GetYaxis()->GetBinLowEdge(j);
	      double norm = (j<3 ? (j==1 ? norm_mip : 7*norm_mip) : 1);
	      if (i<iminH0) norm *= 1./iminH0;
	      h2->SetBinContent(i, j, h3->GetBinContent(ipt, j, i) * norm);
	      h2->SetBinError(i, j, h3->GetBinError(ipt, j, i) * norm);
      } // for j
    } // for i

    // Project out (e)H and E-hadrons
    TH1D *hH = h2->ProjectionX("Hhadrons",1,iminE0);
    TH1D *hE = h2->ProjectionY("Ehadrons",1,iminH0);

    // Calculate mean response over all
    double r0 = (hE->GetMean()+hH->GetMean());
    double s0 = sqrt(pow(hE->GetRMS(),2) + pow(hH->GetRMS(),2));

    TH1D *h = new TH1D("h",";H / T (GeV);E / T (GeV);",75,0,1.5);//100,0,2);
    h->SetMinimum(0);
    h->SetMaximum(1.5);//2);
    
    double ptx = 45;
    // Before photon fix
    double mine = max(0.0,minE);//mip_ecal/pt1;
    double minh = max(0.0,minH);//max(0.1,minpt_hcal/pt1);
    double maxe = min(1.0,1.0-minh);
    double maxh = min(1.0,1.0-mine);
    // Instantiate these only once so we can recycle fit from previous round
    f2 = (f2 ? f2 : new TF2("f2",func,minh,maxh,mine,maxe,6));
    f2b = (f2b ? f2b : new TF2("f2b",func,0.,1,0.,1.,6));
    
    if (!goodfit) 
      f2->SetParameters(3.14, 0.866, 0.0970, 0.638, 0.577, 0.351); // 48-50 GeV
    if (pt1<25) {
      f2->FixParameter(4, 1.); // low pT showers all in ECAL
    }
    
    h2->Fit(f2,ipt0 ? "RN" : "QRN");
    h2->Fit(f2,ipt0 ? "RNME" : "QRNME");
    f2->SetContour(4);
    f2b->SetContour(4);
    for (int i=0; i!=f2->GetNpar(); ++i)
      f2b->SetParameter(i,f2->GetParameter(i));
    f2b->SetLineColor(kBlue+1);
    
    // Fit also E and (e)H-hadrons
    f1 = (f1 ? f1 : new TF1("f1","gaus",0,1));
    // H-hadrons
    f1->SetRange(max(minH, hH->GetMean()-1.5*hH->GetRMS()),
		 min(1.2,  hH->GetMean()+1.5*hH->GetRMS()));
    hH->Fit(f1,"QRN");
    if (f1->GetParError(1)<0.1) {
      int n = ghr->GetN();
      ghr->SetPoint(n, pt, f1->GetParameter(1));
      ghr->SetPointError(n, dpt, f1->GetParError(1));
      gh4->SetPoint(n, 4*pt, f1->GetParameter(1));
      gh4->SetPointError(n, 4*dpt, f1->GetParError(1));
      ghs->SetPoint(n, pt, f1->GetParameter(2) / f1->GetParameter(1));
      ghs->SetPointError(n, dpt, f1->GetParError(2) / f1->GetParameter(1));
    }
    // E-hadrons
    f1->SetRange(max(minE, hE->GetMean()-2*hE->GetRMS()),
		 min(1.,   hE->GetMean()+2*hE->GetRMS()));
    hE->Fit(f1,"QRN");
    if (f1->GetParError(1)<0.1) {
      int n = ger->GetN();
      ger->SetPoint(n, pt, f1->GetParameter(1));
      ger->SetPointError(n, dpt, f1->GetParError(1));
      ges->SetPoint(n, pt, f1->GetParameter(2) / f1->GetParameter(1));
      ges->SetPointError(n, dpt, f1->GetParError(2) / f1->GetParameter(1));
    }
    if (r0>0) {
      int n = gar->GetN();
      gar->SetPoint(n, pt, r0);
      gar->SetPointError(n, dpt, f2->GetParError(1));
      gas->SetPoint(n, pt, s0);
      gas->SetPointError(n, dpt, f2->GetParError(2));
    }

    double rh = f2->GetParameter(1);
    double re = f2->GetParameter(3);
    double ef = f2->GetParameter(4);
    
    if (ipt0==0) {
      // Animated gif scanning from high pT to low pT
      // Remove gif before starting to make sure it won't get messed up?

      int n = grh->GetN();
      if (f2->GetParError(1)<0.1) {
	      goodfit = true;
	      double chi = sqrt(f2->GetChisquare()/f2->GetNDF());
	      double b = f2->GetParameter(1);
	      double a = f2->GetParameter(3);
	      double r = 0.5*(a+b);
	      grh->SetPoint(n, pt, b);
	      grh->SetPointError(n, dpt, f2->GetParError(1)*chi);
	      gre->SetPoint(n, pt, a);
	      gre->SetPointError(n, dpt, f2->GetParError(3)*chi);
	      //
	      grr->SetPoint(n, pt, r);
	      grr->SetPointError(n, dpt, 0.5*(f2->GetParError(1)+
	      				f2->GetParError(3))*chi);
	      grd->SetPoint(n, pt, 0.5*(b-a));
	      grd->SetPointError(n, dpt, 0.5*(f2->GetParError(1)+
	      				f2->GetParError(3))*chi);
	      //
	      grs->SetPoint(n, pt, sqrt(2.)*f2->GetParameter(2) / r );
	      grs->SetPointError(n, dpt, sqrt(2.)*f2->GetParError(2) / r * chi);
	      //
	      gfm->SetPoint(n, pt, f2->GetParameter(4) );
	      gfm->SetPointError(n, dpt, f2->GetParError(4) * chi);
	      gfs->SetPoint(n, pt, f2->GetParameter(4) );
	      gfs->SetPointError(n, dpt, f2->GetParError(5) * chi);
      }

      delete h;
      delete h2;
      delete px;
      delete py;
    }
  } // for ipt

  if (ipt0==0) {

    TH1D *h = new TH1D("h",";T (GeV);Hadron response;",200,0,200);
    h->GetXaxis()->SetRangeUser(ptmin,200);
    h->SetMinimum(0);
    h->SetMaximum(1.1);
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();

    // R: tdrDraw(ghr, "P", kFullCircle, kRed);
    // R: tdrDraw(ger, "P", kFullSquare, kBlue);

    // R: tdrDraw(grh, "P", kOpenTriangleDown, kMagenta+1);//kRed);
    // R: tdrDraw(gre, "P", kOpenTriangleUp, kCyan+1);//kBlue);
    // R: grr->SetMarkerSize(2.0);
    // R: grd->SetMarkerSize(2.0);

    // Use same function as in 1D fit
    TF1 *fhr = new TF1("fhr","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",4,200);
    fhr->SetParameters(1.0806,-0.5,-0.15,-0.3);
    fhr->FixParameter(0, 1.0806); // constrain to old VLE+ME sample
    fhr->SetParLimits(2,-0.18,-0.13);
    cout << "*** Fitting fhr ***" << endl;
    ghr->Fit(fhr,"RN");

    TF1 *frh = new TF1("frh","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    frh->SetParameters(0.9*fhr->GetParameter(0), fhr->GetParameter(1),
		       fhr->GetParameter(2), fhr->GetParameter(3));
    frh->FixParameter(1,fhr->GetParameter(1));
    frh->SetParLimits(2,-0.18,-0.13);
    frh->SetLineColor(kMagenta+1);//Red);
    TF1 *fr = new TF1("fr","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    fr->SetParameters(0.95*frh->GetParameter(0), frh->GetParameter(1),
          frh->GetParameter(2), frh->GetParameter(3));
    fr->SetParLimits(2,-0.18,-0.13);
    fr->SetLineColor(kMagenta+1);
    cout << "*** Fitting frh ***" << endl;
    grh->Fit(frh,"RN");

    TF1 *fre = new TF1("fre","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    fre->SetParameters(0.95*fr->GetParameter(0), fr->GetParameter(1),
		       fr->GetParameter(2), fr->GetParameter(3));
    fre->SetParLimits(2,-0.18,-0.13);
    fre->SetLineColor(kCyan+1);
    cout << "*** Fitting fre ***" << endl;
    gre->Fit(fre,"RN");

    TF1 *fer = new TF1("fer","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",3,200);
    fer->SetParameters(0.98*fre->GetParameter(0), fre->GetParameter(1),
		       fre->GetParameter(2), fre->GetParameter(3));
    fre->SetParLimits(2,-0.18,-0.13);
    fer->SetLineColor(kBlue);
    cout << "*** Fitting fer ***" << endl;
    ger->Fit(fer,"RN");

    cout << endl;
    cout << "New feh: r+/-dr" << endl;
    cout <<Form("feh->SetParameters(0.5, 0.5,\n"
		"                   %1.4f, %1.4f, %1.4f, %1.2f,\n"
		"                   %1.4f, %1.4f, %1.4f, %1.2f);\n",
		frh->GetParameter(0), frh->GetParameter(1),
		frh->GetParameter(2), frh->GetParameter(3),
		fre->GetParameter(0), fre->GetParameter(1),
		fre->GetParameter(2), fre->GetParameter(3));
    cout << endl;
  }

} // draw2D


