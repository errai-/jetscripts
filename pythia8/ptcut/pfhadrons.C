#define pfhadrons_cxx
#include "pfhadrons.h"
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

//#include "tdrstyle_mod12.C";
#include "tdrstyle_mod14.C"
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

inline double oplus(double a, double b) {
  return (sqrt(a*a + b*b));
}

// Double-sided Crystal Ball
Double_t fCrystalBall2(Double_t *xx, Double_t *p) {
  
  double x = xx[0];
  double N0 = p[0]; // Overall normalization
  double alpha1 = p[1]; // start of non-Gaussian tail in units of sigma
  double a1 = fabs(alpha1);
  double n1 = p[2]; // ??
  double alpha2 = p[3];
  double a2 = fabs(alpha2);
  double n2 = p[4];
  double xbar = p[5]; // core Gaussian mean
  double sigma = p[6]; // core Gaussian width

  double A1 = pow(n1 / a1, n1) * exp(-a1*a1/2);
  double B1 = n1 / a1 - a1;
  double C1 = n1 / a1 * 1 / (n1-1) * exp(-a1*a1/2);
  double D1 = sqrt(TMath::Pi()/2) * (1 + TMath::Erf(a1/sqrt(2)));
  double N1 = 1. / (sigma * (C1 + D1));

  double A2 = pow(n2 / a2, n2) * exp(-a2*a2/2);
  double B2 = n2 / a2 - a2;
  double C2 = n2 / a2 * 1 / (n2-1) * exp(-a2*a2/2);
  double D2 = sqrt(TMath::Pi()/2) * (1 + TMath::Erf(a2/sqrt(2)));
  double N2 = 1. / (sigma * (C2 + D2));
		   
  if ( (x-xbar)/sigma > -alpha1 && x-xbar<=0) {
    return ( N0 * N1 * exp(-pow(x-xbar,2)/(2*sigma*sigma)) );
  }
  if ( (x-xbar)/sigma <= -alpha1) {
    return ( N0 * N1 * A1 * pow(B1 - (x-xbar)/sigma, -n1) );
  }
  if ( (x-xbar)/sigma < alpha2 && x-xbar>0) {
    return ( N0 * N2 * exp(-pow(x-xbar,2)/(2*sigma*sigma)) );
  }
  if ( (x-xbar)/sigma >= alpha2) {
    return ( N0 * N2 * A2 * pow(B2 - (xbar-x)/sigma, -n2) );
  }
  
  //assert(false);
  return 0;
}

TH2D *pfhadrons::normHisto(TH2D *h2) {

  TH2D *h2n = (TH2D*)h2->Clone(Form("%sn",h2->GetName()));
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    
    double dx = 1;//h2->GetXaxis()->GetBinWidth(i);
    TH1D *h1 = h2->ProjectionY("_py",i,i);
    double norm = h1->Integral();
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      double dy = h1->GetBinWidth(j);
      if (norm!=0) {
	h2n->SetBinContent(i, j, h2->GetBinContent(i,j) / norm / dx / dy);
	h2n->SetBinError(i, j, h2->GetBinError(i,j) / norm / dx / dy);
      }
      else {
	h2n->SetBinContent(i, j, 0);
	h2n->SetBinError(i, j, 0);
      }
    } // for j
    delete h1;
  } // for i

  return h2n;
} // normHisto

// H+E = f(T; H, E);
// x - variables, p - parameters
Double_t func_eh(Double_t *x, Double_t *p) {

  double T = x[0];
  double E = p[0];
  double H = p[1];
  //Double_t *pa = &p[2];
  //Double_t *pb = &p[5];
  //double a = max(0.33, pa[2]-pa[0]*pow(T,-pa[1]));
  //double b = max(0.33, pb[2]-pb[0]*pow(T,-pb[1]));
  // V2 with r+/-dr
  //Double_t *pr = &p[2];
  //Double_t *pd = &p[6];
  //double r = pr[0]+pr[1]*pow(T,pr[2]);
  //double dr = pd[0]+pd[1]*pow(T,pd[2]);
  //double a = max(0.33,r - dr);
  //double b = max(0.33,r + dr);
  Double_t *pa = &p[6];
  Double_t *pb = &p[2];
  double a = max(0.33, (pa[0]+pa[1]*pow(T+pa[3],pa[2])) * max(0., T+pa[3])/T);
  double b = max(0.33, (pb[0]+pb[1]*pow(T+pb[3],pb[2])) * max(0., T+pb[3])/T);

  double EH = a*(T - H/b) + b*(T - E/a);

  return EH;
} // func_he

void pfhadrons::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L pfhadrons.C
//      Root > pfhadrons t
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

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   int n1 = 32469724;
   int n2 = 3695873;
   int ntot = n1+n2;
   int ntot2 = 1208951;
   int ntot3 = ntot2 + n2;
   assert(nentries==ntot || nentries==ntot2 || nentries==ntot3);
   //nentries = 1000000;// 3%
   //nentries = 5000000;// 15%

   TFile *fout = new TFile("pfhadrons.root","recreate");

   double x[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
		 1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
		 2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
		 3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
		 4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,// 0.1
		 5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,
		 7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,
		 9.0,9.2,9.4,9.6,9.8,
		 10.0,10.2,10.4,10.6,10.8,               // 0.2
		 11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,
		 16,16.5,17,17.5,18,18.5,19,19.5,        // 0.5
		 20,21,22,23,24,25,26,27,28,29,
		 30,31,32,33,34,35,36,37,38,39,
		 40,41,42,43,44,45,46,47,48,49,          // 1.0
		 50,52,54,56,58,60,62,64,66,68,
		 70,72,74,76,78,80,82,84,86,88,
		 90,92,94,96,98,                         // 2.0
		 100,105,110,115,120,125,130,135,140,145,
		 150,155,160,165,170,175,180,185,190,195,
		 200                                     // 5.0
   };
   const int nx = sizeof(x)/sizeof(x[0])-1;

   double xw[] = {0.0,0.2,0.4,0.6,0.8,
		  1.0,1.2,1.4,1.6,1.8,
		  2.0,2.2,2.4,2.6,2.8,
		  3.0,3.2,3.4,3.6,3.8,
		  4.0,4.2,4.4,4.6,4.8,// 0.2
		  5.0,5.4,5.8,6.2,6.6,
		  7.0,7.4,7.8,8.2,8.6,
		  9.0,9.2,9.4,9.6,9.8, // 0.4
		  10,10.4,             // 0.4-0.6
		  11,12,13,14,15,
		  16,17,18,19,        // 1.0
		  20,22,24,26,28,
		  30,32,34,36,38,
		  40,42,44,46,48,     // 2.0
		  50,54,58,62,66,
		  70,74,78,82,86,     // 4.0
		  90,94,              // 4-6
		  100,110,120,130,140,
		  150,160,170,180,190,
		  200                 // 10
   };
   const int nxw = sizeof(xw)/sizeof(xw[0])-1;
		 
   const int nz = 100;
   double z[nz+1];
   for (int i = 0; i != nz+1; ++i) z[i] = 2./nz*i;

   // H-hadrons (ECAL MIPs with no observable deposit)
   TH2D *h2r = new TH2D("h2r","H-hadrons;T (GeV);H/T",nx,x,500,0,5);
   TH2D *h2c = new TH2D("h2c","H-hadrons;T (GeV);Hc/T",nx,x,500,0,5);

   // eH-hadrons (ECAL MIPs)
   TH2D *he2r = new TH2D("he2r","eH-hadrons;T (GeV);(e+H)/T",nx,x,500,0,5);
   TH2D *he2x = new TH2D("he2x","eH-hadrons;T (GeV);(e+Hc)/T",nx,x,500,0,5);
   TH2D *he2c = new TH2D("he2c","eH-hadrons;T (GeV);(e+H)c/T",nx,x,500,0,5);

   // EH-hadrons (shower starting in ECAL)
   TH2D *eh2r = new TH2D("eh2r","EH-hadrons;T (GeV);(E+H)/T",nx,x,500,0,5);
   TH2D *eh2x = new TH2D("eh2x","EH-hadrons;T (GeV);(E+Hc)/T",nx,x,500,0,5);
   TH2D *eh2c = new TH2D("eh2c","EH-hadrons;T (GeV);(E+H)c/T",nx,x,500,0,5);

   // E-hadrons (ECAL only deposits)
   TH2D *e2r = new TH2D("e2r","E-hadrons;T (GeV);E/T",nx,x,500,0,5);
   TH2D *e2c = new TH2D("e2c","E-hadrons;T (GeV);Ec/T",nx,x,500,0,5);

   // All hadrons
   TH2D *a2r = new TH2D("a2r","All hadrons;T (GeV);(E+H)/T",nx,x,500,0,5);
   TH2D *a2x = new TH2D("a2x","All hadrons;T (GeV);(E+Hc)/T",nx,x,500,0,5);
   TH2D *a2c = new TH2D("a2c","All hadrons;T (GeV);(E+H)c/T",nx,x,500,0,5);
   TH3D *e3r = new TH3D("e3r","All hadrons;T (GeV);E/T;H/T",nxw,xw,nz,z,nz,z);

   // response vs E/(E+H) vs pT for tracking optimal mip threshold
   TProfile2D *p2r = new TProfile2D("p2r",";T (GeV);E;(E+H)/T",nx,x,nx,x);

   TH2D *h2ef = new TH2D("h2ef","ECAL fraction;T (GeV);E/(E+H)",nx,x,100,0,1);

   TProfile *pefa = new TProfile("pefa","ECAL fraction;T (GeV);E/(E+H)",nx,x);
   TProfile *pefb = new TProfile("pefb","ECAL fraction;T (GeV);E/(E+H)",nx,x);
   TProfile *pefh = new TProfile("pefh","ECAL fraction;T (GeV);E/(E+H)",nx,x);
   TProfile *pefe = new TProfile("pefe","ECAL fraction;T (GeV);E/(E+H)",nx,x);
   TProfile *pef0 = new TProfile("pef0","ECAL fraction;T (GeV);E/(E+H)",nx,x);
   TProfile *pef1 = new TProfile("pef1","ECAL fraction;T (GeV);E/(E+H)",nx,x);

   TProfile *pfh = new TProfile("pfh","H-hadrons;T (GeV);f(E=0,H>0)",nx,x);
   TProfile *pfe = new TProfile("pfe","E-hadrons;T (GeV);f(E>0,H=0)",nx,x);
   TProfile *pf0 = new TProfile("pf0","0-hadrons;T (GeV);f(E=0,H=0)",nx,x);
   TProfile *pfeh = new TProfile("pfeh","EH-hadrons;T (GeV);f(E>0,H>0)",nx,x);
   TProfile *pfe0h = new TProfile("pfe0h","EH-hadrons;T (GeV);f(E<1,H>0)",nx,x);
   TProfile *pfe1h = new TProfile("pfe1h","EH-hadrons;T (GeV);f(E>1,H>0)",nx,x);

   //TF1 *fp = new TF1("fp","x*([0]+[1]*pow(max(x,[4]),[2])+[3]/max(x,[4]))",
   //	     0,4000);
   // R(T) constant below point where R(T)*T=0.8 GeV => T<3.7 GeV
   //fp->SetParameters(1.25191,-0.296681, -0.0598,-2.87122,3.7);

   // New fit function taking into account energy loss before HCAL
   //TF1 *fp1 = new TF1("fp1","x * ([0]+[1]*pow(max(x,[4]),[2]))"
   // first version was missing +[3] from powerlaw
   TF1 *fp1 = new TF1("fp1","x * ([0]+[1]*pow(max(x,[4])+[3],[2]))"
   		      " * (max(x,[4])+[3])/max(x,[4])",
   		      0,4000);
   //fp1->SetParameters( 1.08063, -0.129273, -0.130000,  -2.40708, 3.23); // H
   //fp1->SetParameters(1.0806, -0.0844, -0.1300, -2.4212, 3.62); // H
   fp1->SetParameters(1.0806, -0.0838, -0.1300, -2.3051, 3.44); // He

   // After-thought: we should limit true response to 0.33, because at low pT
   // (pT<6 GeV?) we have at most one full nuclear interaction, of which
   // 1/3 goes to pi0 with R=1, and 2/3 goes to pi+// with R>=0

   // "After-burner" to account for net bias from JER+non-linear correction
   // (seems small after fp->fp1, but do anyway)
   //TF1 *fp2 = new TF1("fp2","x * ([0]+[1]*pow(max(x,[4]),[2]))"
   //		      " * (max(x,[4])+[3])/max(x,[4])",
   //		      0,4000);
   //fp2->SetParameters(1.00890, -0.016691, -0.13, 0.0337466, 3.23);

   // EH-hadrons
   TF1 *feh = new TF1("feh",func_eh,0,4000,10);
   //feh->SetParameters(0.5,0.5,
   //		      1.1116, 0.2500, 1.1000,
   //		      0.8195, 0.2499, 1.1000);
   // V2 with r+/-dr
   //feh->SetParameters(0.5, 0.5,
   //	      1.3059, -1.0486, -0.1574,
   //	      -0.0561, 0.2279, -0.1830);
   // V3 with rh, re
   feh->SetParameters(0.5, 0.5,
		      1.4160, -1.0004, -0.1800, -2.19,
		      0.9656, -0.5652, -0.1800, -0.18);

   // E-hadrons
   TF1 *fe = new TF1("fe","x * ([0]+[1]*pow(max(x,[4])+[3],[2]))"
		     " * (max(x,[4])+[3])/max(x,[4])",
		     0,4000);
   fe->SetParameters(1.5704, -1.4278, -0.1004, 0.0000, 4.10);

   Long64_t nbytes = 0, nb = 0; int cnt(0);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t kentry = jentry;
     if (nentries<ntot && nentries!=ntot2 && nentries!=ntot3)
       kentry = (jentry<nentries/2 ? jentry : jentry + n1);

      Long64_t ientry = LoadTree(kentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(kentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (++cnt%100000==0) cout << "." << flush;

      // H-hadrons
      if (fabs(eta)<1.3 && ecal==0) {
	h2r->Fill(p, hcal/p);
	// zero maps to zero so keep lower range 0
	double hc = fp1->GetX(hcal, 0., 4000);
	h2c->Fill(p, hc/p);
      }
      // eH-hadrons (including H and e as subsets)
      if (fabs(eta)<1.3 && ecal<mip_ecal) {
	he2r->Fill(p, (hcal+ecal)/p);
	double hc = fp1->GetX(hcal, 0., 4000);
	he2x->Fill(p, (hc+ecal)/p);
	double hec = fp1->GetX(hcal+ecal, 0., 4000);
	he2c->Fill(p, hec/p);
      }
      // EH-hadrons
      if (fabs(eta)<1.3 && ecal>=mip_ecal && hcal>0) {
	eh2r->Fill(p, (hcal+ecal)/p);
	double ehx = fp1->GetX(hcal+ecal, 0., 4000.);
	eh2x->Fill(p, ehx/p);
	feh->SetParameter(0,ecal);
	feh->SetParameter(1,hcal);
	double ehc = feh->GetX(hcal+ecal, 0., 4000.);
	eh2c->Fill(p, ehc/p);
      }
      // E-hadrons not dealt with yet
      if (fabs(eta)<1.3 && ecal>mip_ecal && hcal==0) {
	e2r->Fill(p, ecal/p);
	double ec = fe->GetX(ecal, 0., 4000.);
	e2c->Fill(p, ec/p);
      }
      // All hadrons
      if (fabs(eta)<1.3) {
	
	a2r->Fill(p, (ecal+hcal)/p);
	double hc = fp1->GetX(hcal, 0., 4000.);
	a2x->Fill(p, (ecal+hc)/p);

	if (ecal<mip_ecal) { // H, eH and lost hadrons
	  double hec = fp1->GetX(hcal+ecal, 0., 4000);
	  a2c->Fill(p, hec/p);
	}
	if (ecal>=mip_ecal && hcal>0) { // EH-hadrons
	  feh->SetParameter(0,ecal);
	  feh->SetParameter(1,hcal);
	  double ehc = feh->GetX(hcal+ecal, 0., 4000.);
	  a2c->Fill(p, ehc/p);
	}
	if (ecal>mip_ecal && hcal==0) { // E-hadrons
	  double ec = fe->GetX(ecal, 0., 4000.);
	  a2c->Fill(p, ec/p);
	}

	if (ecal+hcal>0)
	  p2r->Fill(p, ecal, (ecal+hcal)/p);
	
	// Ultimate correlation plot of ECAL and HCAL responses
	e3r->Fill(p, ecal/p, hcal/p);
      }


      // ECAL fractions
      if (fabs(eta)<1.3 && ecal>0 && hcal>0) {
	pefb->Fill(p, ecal/(ecal+hcal));
	if (ecal<mip_ecal) pef0->Fill(p, ecal/(ecal+hcal)); // eH
	if (ecal>mip_ecal) pef1->Fill(p, ecal/(ecal+hcal)); // H
      }
      if (fabs(eta)<1.3 && hcal>0) {
	pefh->Fill(p, ecal/(ecal+hcal));
      }
      if (fabs(eta)<1.3 && ecal>0) {
	pefe->Fill(p, ecal/(ecal+hcal));
      }
      if (fabs(eta)<1.3 && (ecal>0 || hcal>0)) {
	pefa->Fill(p, ecal/(ecal+hcal));
	h2ef->Fill(p, ecal/(ecal+hcal));
      }
      if (fabs(eta)<1.3) {
	pfh->Fill(p, (ecal==0 && hcal>0 ) ? 1 : 0);
	pfe->Fill(p, (ecal>0  && hcal==0) ? 1 : 0);
	pf0->Fill(p, (ecal==0 && hcal==0) ? 1 : 0);
	pfeh->Fill(p, (ecal>0 && hcal>0)  ? 1 : 0);
	pfe1h->Fill(p, (ecal>1 && hcal>0)  ? 1 : 0);
	pfe0h->Fill(p, (ecal>0 && ecal<mip_ecal && hcal>0)  ? 1 : 0);
      }
      /*
      if (fabs(eta)<1.3 && ecal!=0 && hcal!=0) {
	double hc1 = fp1->GetX(hcal, 0., 4000);
	double hc2 = fp2->GetX(hc1, 0., 4000);
	eh2r->Fill(p, (hc2+ecal)/p);
	eh2h->Fill(p, hc2/p);
      }
      */
   }

   h2r->Write();
   h2c->Write();

   he2r->Write();
   he2x->Write();
   he2c->Write();

   eh2r->Write();
   eh2x->Write();
   eh2c->Write();

   e2r->Write();
   e2c->Write();

   a2r->Write();
   a2c->Write();
   e3r->Write();

   p2r->Write();
   h2ef->Write();

   pefa->Write();
   pefb->Write();
   pefh->Write();
   pefe->Write();
   pef0->Write();
   pef1->Write();

   pfh->Write();
   pfe->Write();
   pf0->Write();
   pfeh->Write();
   pfe0h->Write();
   pfe1h->Write();

   normHisto(h2r)->Write();
   normHisto(h2c)->Write();

   normHisto(eh2r)->Write();
   normHisto(eh2x)->Write();
   normHisto(eh2c)->Write();

   normHisto(he2r)->Write();
   normHisto(he2x)->Write();
   normHisto(he2c)->Write();

   normHisto(e2r)->Write();
   normHisto(e2c)->Write();

   normHisto(a2r)->Write();
   normHisto(a2c)->Write();

   normHisto(h2ef)->Write();
   
   fout->Close();
} // Loop


void pfhadrons::draw1D(string mode) {
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();
  gStyle->SetPalette(1);

  TFile *f = new TFile("pfhadrons.root","read");
  assert(f && !f->IsZombie());
  
  double ptmin = minpt_hcal;
  double ptmax = 200;
  TH2D *h2(0);
  if (mode=="H")    { h2 = (TH2D*)f->Get("h2rn"); }
  if (mode=="Hc")   { h2 = (TH2D*)f->Get("h2cn");  ptmin = minpt_hcal/0.33; }
  if (mode=="He")   { h2 = (TH2D*)f->Get("he2rn"); ptmin = minpt_ecal; }
  if (mode=="Hex")  { h2 = (TH2D*)f->Get("he2xn"); ptmin = minpt_ecal; }
  if (mode=="Hec")  { h2 = (TH2D*)f->Get("he2cn"); ptmin = minpt_ecal/0.33; }
  if (mode=="EH")   { h2 = (TH2D*)f->Get("eh2rn"); ptmin = minpt_ecal; }
  if (mode=="EHx")  { h2 = (TH2D*)f->Get("eh2xn"); ptmin = minpt_ecal/0.33; }
  if (mode=="EHc")  { h2 = (TH2D*)f->Get("eh2cn"); ptmin = minpt_ecal/0.33; }
  if (mode=="E")    { h2 = (TH2D*)f->Get("e2rn");  ptmin = minpt_ecal; }
  if (mode=="Ec")   { h2 = (TH2D*)f->Get("e2cn");  ptmin = minpt_ecal/0.33; }
  if (mode=="A")    { h2 = (TH2D*)f->Get("a2rn");  ptmin = minpt_ecal; }
  if (mode=="Ac")   { h2 = (TH2D*)f->Get("a2cn");  ptmin = minpt_ecal/0.33; }
  assert(h2);

  // Clean out lost hadrons, single spikes
  if (mode=="Ec") {
    for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	if (h2->GetXaxis()->GetBinCenter(i)<1.5) {
	  h2->SetBinContent(i, j, 0);
	  h2->SetBinError(i, j, 0);
	}
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetRightMargin(0.15);
  c1->SetLogx();
  c1->SetLogz();
  //h2->GetXaxis()->SetRangeUser(ptmin,200);
  h2->GetXaxis()->SetRangeUser(0.1,200);
  h2->GetYaxis()->SetRangeUser(0.01,3.01); // 300 bins
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->Draw("COLZ");
  
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  //TH1D *h = new TH1D("h","H-hadrons;H/T;Normalized area",500,0,5);
  //h->DrawClone("AXIS");

  // mean
  TGraphErrors *glmean = new TGraphErrors(0);
  TGraphErrors *glmode = new TGraphErrors(0);
  TGraphErrors *glmed = new TGraphErrors(0);
  TGraphErrors *g1 = new TGraphErrors(0);
  TGraphErrors *gg = new TGraphErrors(0);
  TGraphErrors *gc = new TGraphErrors(0);
  TGraphErrors *ga = new TGraphErrors(0);

  // resolution
  TGraphErrors *ggs = new TGraphErrors(0);
  TGraphErrors *gcs = new TGraphErrors(0);

  // Log-normal function as described in Wikipedia
  TF1 *fl = new TF1(Form("fl_%d",0),"[2]*1/(x*[1]*sqrt(2*TMath::Pi()))"
		    "*exp(-pow(log(x)-[0],2)/(2*[1]*[1]))",0,5);
  fl->SetParameters(-0.56, 0.62, 0.67);
  

  // ROOT's version of LogNormal with extra shape parameter
  TF1 *f1 = new TF1(Form("f1_%d",0),
		    "[3]*TMath::LogNormal(x,[0],[1],[2])",0,5);
  f1->SetParameters(0.708, 0.023, 0.532, 1);
  //f1->SetParameters(-0.56, 0, 0.62, 0.67);

  TF1 *fc = new TF1(Form("fc_%d",0),fCrystalBall2,0,5,7);
  fc->SetParameters(1, 100,0, 1,2, 1,0.1); // 1-sided
  fc->FixParameter(1, 100); fc->FixParameter(2, 0); // 1-sided
  //fc->SetParameters(1, 1,2, 1,2, 1,0.1); // 2-sidd
  //fc->SetParLimits(1, 1,10);
  fc->SetParLimits(3, 1,10);
  bool isgoodfitc = false;

  //for (int ipt = 10; ipt != 2000; ++ipt) {
  bool loop = true;
  
  int iptn = h2->GetXaxis()->FindBin(ptmax-1e-3);
  for (int ipt = iptn; ipt != 0; --ipt) {
    //{  int ipt = 490;

    //int ipt = 100; 
    int dipt = 1;
    int jpt = ipt + dipt;
    double pt1 = h2->GetXaxis()->GetBinLowEdge(ipt);
    double pt2 = h2->GetXaxis()->GetBinLowEdge(jpt);
    double dpt = (pt2-pt1);

    TH1D *h1 = h2->ProjectionY("_py",ipt,jpt-1);
    h1->Scale(1./dipt);
    //h1->Rebin(4); h1->Scale(1./4.);
    TH1D *h0 = (TH1D*)h1->Clone("h0");
    h0->SetBinContent(1,0); h0->SetBinError(1,0);
    if (!loop) h0->DrawClone();//"SAME");
    
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.045);
    tex->SetNDC();
    tex->DrawLatex(0.45,0.87,Form("%1.1f < T < %1.1f GeV",pt1,pt2));
    
    double xmin = h1->GetBinLowEdge(h1->FindBin(ptmin/pt1)+1);

    fl->SetRange(xmin,5);
    h1->Fit(fl,loop ? "QRN" : "RN");
    fl->SetLineStyle(kSolid);
    if (!loop) fl->DrawClone("SAME");
    fl->SetRange(0,5);
    fl->SetLineStyle(kDashed);
    if (!loop) fl->DrawClone("SAME");
    
    tex->SetTextColor(kRed);
    tex->DrawLatex(0.45,0.73,Form("#chi^{2} / NDF = %1.1f / %d",
				  fl->GetChisquare(), fl->GetNDF()));
    
    f1->SetParameters(fl->GetParameter(1), 0, exp(fl->GetParameter(0)),
		      fl->GetParameter(2));
    //f1->SetParameters(0.708, 0.023, 0.532, 1);
    f1->SetLineStyle(kSolid);
    f1->SetLineColor(kBlue);
    f1->SetRange(xmin,5);
    h1->Fit(f1,loop ? "QRN" : "RN");
    if (!loop) f1->DrawClone("SAME");
    f1->SetRange(0,5);
    f1->SetLineStyle(kDashed);
    if (!loop) f1->DrawClone("SAME");
    
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.45,0.80,Form("#chi^{2} / NDF = %1.1f / %d",
				  f1->GetChisquare(), f1->GetNDF()));
    
    // Regular Gaussian fit
    TF1 *fg = new TF1(Form("fg_%d",ipt),"gaus",xmin,5);
    h1->Fit(fg,loop ? "QRN" : "RN");
    fg->SetLineColor(kGreen+2);
    if (!loop) fg->DrawClone("SAME");
    fg->SetRange(0,5);
    fg->SetLineStyle(kDashed);
    if (!loop) fg->DrawClone("SAME");
    
    tex->SetTextColor(kGreen+1);
    tex->DrawLatex(0.45,0.66,Form("#chi^{2} / NDF = %1.1f / %d",
				  fg->GetChisquare(), fg->GetNDF()));
    
    // Crystal Baal fit
    fc->SetRange(xmin,5);
    //fc->SetParameters(fg->GetParameter(0), 100,0,
    //	      fc->GetParameter(3), fc->GetParameter(4),
    //	      fg->GetParameter(1), fg->GetParameter(2));
    //fc->SetParLimits(1, 1,5);
    if (!isgoodfitc) {
      fc->SetParameters(fg->GetParameter(0), 100,0,
			1, 5,
			fg->GetParameter(1), fg->GetParameter(2));
    }
    //fc->SetParLimits(1, 1,5);
    //fc->SetParLimits(3, 1,5);
    h1->Fit(fc,loop ? "QRN" : "RN");
    fc->SetLineStyle(kSolid);
    fc->SetLineColor(kBlue);
    if (!loop) fc->DrawClone("SAME");
    fc->SetRange(0,5);
    fc->SetLineStyle(kDashed);
    if (!loop) fc->DrawClone("SAME");

    h1->Rebin(10); h1->Scale(1./10.);
    h1->SetLineColor(kBlack);
    h1->DrawClone("SAME HIST");
    
    //c2->SaveAs(Form("pdf/pfhadrons_fit%1.1f-%1.1f.pdf",0.1*ipt,0.1*jpt));

    double pt = 0.5*(pt1+pt2);//0.1*0.5*(ipt+jpt);
    double ept = 0.5*(pt2-pt1);//0.1*0.5*(jpt-ipt);

    if (fg->GetParError(1)<0.1 && h1->Integral()>0) {
      int n = gg->GetN();
      gg->SetPoint(n, pt, fg->GetParameter(1));
      gg->SetPointError(n, ept, fg->GetParError(1));
      ggs->SetPoint(n, pt, fg->GetParameter(2));
      ggs->SetPointError(n, ept, fg->GetParError(2));
    }

    if (fc->GetParError(5)<0.1 && h1->Integral()>0) {
      isgoodfitc = true;
      int n = gc->GetN();
      gc->SetPoint(n, pt, fc->GetParameter(5));
      gc->SetPointError(n, ept, fc->GetParError(5));
      gcs->SetPoint(n, pt, fc->GetParameter(6));
      gcs->SetPointError(n, ept, fc->GetParError(6));
    }
    else {
      isgoodfitc = false;
    }

    if (h1->Integral()>0) {
      int n = ga->GetN();
      ga->SetPoint(n, pt, h1->GetMean());
      ga->SetPointError(n, ept, h1->GetMeanError());
    }

    // http://en.wikipedia.org/wiki/Log-normal_distribution
    // mu = log(m*m/sqrt(v+m*m)),  sigma = sqrt(log(1+v/(m*m)));
    // => m = E[X] = exp(mu + 0.5*sigma*sigma)
    // => v = SD[x] = m*sqrt(exp(sigma*sigma) - 1);
    // Mode[X] = exp(mu-sigma*sigma); Med[X] = exp(mu);
    if (fl->GetParError(0)<0.1) {
      int n = glmean->GetN();
      double mu = fl->GetParameter(0);
      double sigma = fl->GetParameter(1);
      double mean = exp(mu + 0.5*sigma*sigma);
      double mode = exp(mu - sigma*sigma);
      double med = exp(mu);
      glmean->SetPoint(n, pt, mean);
      glmean->SetPointError(n, ept, 0);//fg->GetParError(1)); // fix
      glmode->SetPoint(n, pt, mode);
      glmode->SetPointError(n, ept, 0);//fg->GetParError(1)); // fix
      glmed->SetPoint(n, pt, med);
      glmed->SetPointError(n, ept, 0);//fg->GetParError(1)); // fix
    }
  } // for ipt

  c1->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0.1,1,200,1);

  ga->SetLineColor(kBlack);
  ga->DrawClone("SAMEL");

  gc->SetLineColor(kBlue);
  gc->DrawClone("SAMEL");

  gg->SetLineColor(kGreen+2);
  gg->DrawClone("SAMEL");

  glmean->SetLineStyle(kDashed);
  glmean->SetLineColor(kRed);
  glmean->Draw("SAMEL");
  glmode->SetLineColor(kRed);
  glmode->SetLineStyle(kDotted);
  glmode->Draw("SAMEL");
  glmed->SetLineColor(kRed);
  glmed->SetLineStyle(kSolid);
  glmed->Draw("SAMEL");

  // for fitting, add small "shape uncertainty" to graph
  //for (int i = 0; i != gg->GetN(); ++i) {
  //gg->SetPointError(i, gg->GetEX()[i], oplus(gg->GetEY()[i],0.02));
  //}

  TF1 *fp = new TF1("fp","[0]+[1]*pow(x,[2])+[3]/x",9,200);
  //fp->SetParameters(1,-0.5,-0.5,0);
  fp->SetParameters(1,-0.5,-0.15,-0.3);
  fp->SetParLimits(2,-0.18,-0.13);
  // [2]=m-1, m=0.82-0.87 => [2] ~ -0.15 (-0.18 to -0.13
  // -[1]=(1-h/e)/E0^(m-1), h/e~0.5? E0~1 GeV? => [1] ~ -0.5
  // MIP loss in ECAL ~ 1.5 MeV cm^2/g * 10 g/cm^2 * 22 cm = 0.3 GeV,
  //fp->SetParameters(1.25191,-0.296681, -0.0598,-2.87122);
  gg->Fit(fp,"RNME");
  fp->SetLineColor(kBlue);
  fp->DrawClone("SAME");
  fp->SetLineStyle(kDashed);
  fp->SetRange(1,200);
  fp->DrawClone("SAME");

  cout << Form("fp(50) = %1.3f", fp->Eval(50.)) << endl;
  double dh = fp->GetParameter(3);
  cout << Form("[fp(50') = %1.3f]", fp->Eval(50)-dh/50 ) << endl;
  fp->SetParameter(3,0);
  cout << Form("[fp(50') = %1.3f]", fp->Eval(50) ) << endl;
  cout << Form("[fp(1') = %1.3f]", fp->Eval(1) ) << endl;
  cout << Form("fp chi2/ndf = %1.1f / %d\n",fp->GetChisquare(),fp->GetNDF());

  // Power law response with T shifted by dT (energy loss in ECAL)
  TF1 *fd = new TF1("fd","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
  //fd->SetParameters(1,-0.5,-0.5,0);
  if (mode=="E") {
    fd->SetRange(5,40);
    fd->SetParameters(1.,-1,-0.15,0);
    fd->SetParLimits(2,-2,-0.1);
    fd->FixParameter(3,0);
    gg->Fit(fd,"RNW");
  }
  else if (mode=="He" || mode=="H") {
    fd->SetRange(7,49);
    fd->SetParameters(1.08,-0.13,-0.13,-2.4);
    fd->SetParLimits(2,-0.18,-0.13);
    //fp1->SetParameters( 1.08063, -0.129273, -0.130000,  -2.40708, 3.23);
    fd->FixParameter(0, 1.08063); // Use old results to constrain high pT
    gg->Fit(fd,"RNME");
  }
  else {
    fd->SetParameters(1,-0.5,-0.15,-0.3);
    fd->SetParLimits(2,-0.18,-0.13);
    //fd->SetParameters(1.05,-1,-0.7,-1);
    gg->Fit(fd,"RNME");
  }
  fd->SetLineColor(kMagenta+2);
  fd->DrawClone("SAME");
  fd->SetLineStyle(kDashed);
  fd->SetRange(-fd->GetParameter(3),200);
  fd->DrawClone("SAME");

  cout << Form("fp1->SetParameters(%1.4f, %1.4f, %1.4f, %1.4f, %1.2f);",
	       fd->GetParameter(0), fd->GetParameter(1),
	       fd->GetParameter(2), fd->GetParameter(3),
	       fd->GetX(0.33, fabs(fd->GetParameter(3))+0.5, 200)
	       ) << endl;

  cout << Form("fd(50) = %1.3f", fd->Eval(50.)) << endl;
  double dt = fd->GetParameter(3);
  cout << Form("[fd(50') = %1.3f]", fd->Eval(50+dt)*50/(50+dt) ) << endl;
  fd->SetParameter(3,0);
  cout << Form("[fd(50') = %1.3f]", fd->Eval(50) ) << endl;
  cout << Form("[fd(1') = %1.3f]", fd->Eval(1) ) << endl;
  cout << Form("fd chi2/ndf = %1.1f / %d\n",fd->GetChisquare(),fd->GetNDF());

  c1->RedrawAxis();
  c1->Update();
  c1->SaveAs(Form("pdf/pfhadrons_%s_vsT.pdf",mode.c_str()));

  // draw resolution
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  c3->SetLogx();

  TH1D *h = new TH1D("h",Form(";T (GeV);RMS(%s)",h2->GetYaxis()->GetTitle()),
		     200,1,200);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(1.0);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  ggs->Draw("SAMEP");

  TFile *fout = new TFile("pfhadrons_jer.root","UPDATE");
  ggs->SetMarkerStyle(kOpenCircle);
  ggs->Write(Form("jerg_%s",mode.c_str()),TObject::kOverwrite);
  gcs->SetMarkerStyle(kFullCircle);
  gcs->Write(Form("jerc_%s",mode.c_str()),TObject::kOverwrite);
  //
  gg->SetMarkerStyle(kOpenCircle);
  gg->SetMarkerColor(gg->GetLineColor());
  gg->Write(Form("jesg_%s",mode.c_str()),TObject::kOverwrite);
  gc->SetMarkerStyle(kFullCircle);
  gc->SetMarkerColor(gc->GetLineColor());
  gc->Write(Form("jesc_%s",mode.c_str()),TObject::kOverwrite);
  fout->Close();

} // draw1D

void pfhadrons::drawFracs() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile("pfhadrons.root","READ");
  assert(f && !f->IsZombie());

  TProfile *pf0 = (TProfile*)f->Get("pf0"); assert(pf0);
  TProfile *pfe = (TProfile*)f->Get("pfe"); assert(pfe);
  TProfile *pfh = (TProfile*)f->Get("pfh"); assert(pfh);
  TProfile *pfeh = (TProfile*)f->Get("pfeh"); assert(pfeh);
  TProfile *pfe0h = (TProfile*)f->Get("pfe0h"); assert(pfe0h);
  TProfile *pfe1h = (TProfile*)f->Get("pfe1h"); assert(pfe1h);

  curdir->cd();

  TH1D *hf0 = pf0->ProjectionX();
  TH1D *hfe = pfe->ProjectionX();
  TH1D *hfh = pfh->ProjectionX();
  TH1D *hfeh = pfeh->ProjectionX();
  TH1D *hfe0h = pfe0h->ProjectionX();
  TH1D *hfe1h = pfe1h->ProjectionX();

  TH1D *h = new TH1D("h",";True hadron energy (GeV);Hadron fractions",
		     2000,0,200);
  //1990,1,200);
  h->SetMaximum(1.0);
  h->SetMinimum(0.0);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  // Stack
  //hfe->Add(hf0);
  //hfh->Add(hfe);
  //hfeh->Add(hfh);
  hfe0h->Add(hfe1h);
  hfh->Add(hfe0h);
  //hfh->Add(hfeh);
  hfe->Add(hfh);
  hf0->Add(hfe);

  extraText = "Simulation";
  lumi_8TeV = "Pythia 6 Z2*";
  TCanvas *c1 = tdrCanvas("c1", h, 2, 11, true);
  gPad->SetLogx();

  c1->cd();
  //tdrDraw(hfeh,"H",0,0,1,kMagenta+1,1001,kGreen+1);
  //tdrDraw(hfh,"H",0,0,1,kRed+1,1001,kRed);
  //tdrDraw(hfe,"H",0,0,1,kBlue+1,1001,kBlue);
  //tdrDraw(hf0,"H",0,0,1,kGray+1,1001,kGray);
  tdrDraw(hf0,"H",0,0,1,kGray+1,1001,kGray);
  tdrDraw(hfe,"H",0,0,1,kBlue+1,1001,kBlue);
  tdrDraw(hfh,"H",0,0,1,kRed+1,1001,kRed);
  //tdrDraw(hfeh,"H",0,0,1,kGreen+2,1001,kGreen+1);
  tdrDraw(hfe0h,"H",0,0,1,kMagenta+2,1001,kMagenta+1);
  tdrDraw(hfe1h,"H",0,0,1,kGreen+2,1001,kGreen+1);
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  // 0-50
  /*
  tex->SetTextColor(kGreen+3);
  tex->DrawLatex(0.73,0.41,"EH-hadrons");
  tex->SetTextColor(kRed+3);
  tex->DrawLatex(0.77,0.79,"H-hadrons");
  tex->DrawLatex(0.645,0.45,"H");
  tex->SetTextColor(kBlue+3);
  tex->DrawLatex(0.47,0.50,"E-hadrons");
  tex->SetTextColor(kGray+3);
  tex->DrawLatex(0.19,0.56,"Lost hadrons");
  */
  // 0-200
  //tex->SetTextSize(0.045);
  /*
  tex->SetTextColor(kGray+3);
  tex->DrawLatex(0.18,0.56,"Lost hadrons");
  tex->SetTextColor(kBlue+3);
  tex->DrawLatex(0.42,0.51,"E-had.");
  tex->SetTextColor(kRed+3);
  tex->DrawLatex(0.55,0.46,"H");
  tex->DrawLatex(0.695,0.82,"H-hadrons");
  tex->SetTextColor(kGreen+3);
  tex->DrawLatex(0.68,0.41,"EH-hadrons");
  */
  // 0-200, with EH->eH+EH
  //tex->SetTextSize(0.045);
  tex->SetTextColor(kGray+3);
  tex->DrawLatex(0.18,0.56,"Lost hadrons");
  tex->SetTextColor(kBlue+3);
  tex->DrawLatex(0.415,0.495,"E-had.");
  tex->SetTextColor(kRed+3);
  tex->DrawLatex(0.535,0.43,"H");
  tex->DrawLatex(0.695,0.82,"H-hadrons");
  tex->SetTextColor(kMagenta+3);
  //tex->DrawLatex(0.565,0.365,"eH");
  tex->DrawLatex(0.58,0.38,"eH");
  tex->DrawLatex(0.68,0.56,"eH-hadrons");
  tex->SetTextColor(kGreen+3);
  tex->DrawLatex(0.68,0.30,"EH-hadrons");

  CMS_lumi(c1, 2, 11);
  gPad->RedrawAxis();

  c1->SaveAs("pdf/pfhadrons_frac.pdf");

  TCanvas *c2 = tdrCanvas("c2",h,2,33,true);
  gPad->SetLogx();

  tdrDraw(pf0,"H",0,0,1,kGray+3,0,0);
  tdrDraw(pfe,"H",0,0,1,kBlue+3,0,0);
  tdrDraw(pfh,"H",0,0,1,kRed+3,0,0);
  tdrDraw(pfe0h,"H",0,0,1,kMagenta+3,0,0);
  tdrDraw(pfe1h,"H",0,0,1,kGreen+3,0,0);

  TLegend *leg = tdrLeg(0.60,0.50,0.75,0.75);
  leg->AddEntry(pf0,"Lost");
  leg->AddEntry(pfe,"E");
  leg->AddEntry(pfh,"H");
  leg->AddEntry(pfe0h,"eH");
  leg->AddEntry(pfe1h,"EH");

  c2->SaveAs("pdf/pfhadrons_fracs.pdf");

} // drawFracs


void pfhadrons::drawECAL() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile("pfhadrons.root","READ");
  assert(f && !f->IsZombie());

  TProfile *pefa = (TProfile*)f->Get("pefa"); assert(pefa);
  TProfile *pefb = (TProfile*)f->Get("pefb"); assert(pefb);
  TProfile *pefe = (TProfile*)f->Get("pefe"); assert(pefe);
  TProfile *pefh = (TProfile*)f->Get("pefh"); assert(pefh);
  TProfile *pef1 = (TProfile*)f->Get("pef1"); assert(pef1);

  curdir->cd();

  TH1D *hefa = pefa->ProjectionX();
  TH1D *hefb = pefb->ProjectionX();
  TH1D *hefe = pefe->ProjectionX();
  TH1D *hefh = pefh->ProjectionX();
  TH1D *hef1 = pef1->ProjectionX();

  TH1D *h = new TH1D("h",";True hadron energy (GeV);ECAL fraction",
		    1990,1,200);
  h->SetMaximum(1.0);
  h->SetMinimum(0.0);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  extraText = "Simulation";
  lumi_8TeV = "Pythia 6 Z2*";
  TCanvas *c1 = tdrCanvas("c1", h, 2, 11, true);
  gPad->SetLogx();

  c1->cd();
  //tdrDraw(hef0,"H",0,0,1,kBlue-8,1001,kBlue-9);
  //tdrDraw(hef,"H",0,0,1,kBlue+1,3001,kBlue);

  //tdrDraw(hefa,"H",0,0,1,kBlue+1-9,1001,kBlue-9);
  tdrDraw(hef1,"H",0,0,1,kBlue-9+2,1001,kBlue-9);
  tdrDraw(hefa,"L",0,0,1,kBlue+1);
  // //tdrDraw(hefb,"H",0,0,1,kGreen+1,0,0);
  //tdrDraw(hef0,"H",0,0,1,kMagenta+1,0,0);
  // tdrDraw(hef1,"H",0,0,1,kGreen+1,0,0);
  //tdrDraw(hefe,"L",0,0,1,kBlue+1,0,0);
  //tdrDraw(hefh,"L",0,0,1,kRed+1,0,0);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.60,0.85,"|#eta|<0.3");

  TLegend *leg = tdrLeg(0.60,0.70,0.85,0.80);
  leg->AddEntry(hefa,"All hadrons","L");
  //leg->AddEntry(hefe,"E+EH+eH");
  // //leg->AddEntry(hefb,"EH");
  leg->AddEntry(hef1,"EH only","FL");
  //leg->AddEntry(hef0,"eH");
  //leg->AddEntry(hefh,"H+EH+eH");
  gPad->RedrawAxis();

  c1->SaveAs("pdf/pfhadrons_ECAL.pdf");
} // drawECAL


Double_t func(Double_t *x, Double_t *p) {
  
  double ht = x[0];
  double et = x[1];
  double n = p[0];
  double rh = p[1];
  double sh = p[2];
  double re = p[3];
  //double se = p[4];
  //double rf = p[5];
  //double sf = p[6];
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

  //if (pow(dh/sh,2) + pow(de/se,2) + pow(df/sf,2) > 4)
  //return 0;
   
  return (n*TMath::Gaus(dh,0,sh,kTRUE)*TMath::Gaus(de,0,se,kTRUE)
	  *TMath::Gaus(df,0,sf,kTRUE));
}

void pfhadrons::draw2D(double pt) {

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


  //TF1 *feh = new TF1("feh",func_eh,0,4000,8);
  //feh->SetParameters(0.5, 0.5,
  //	     1.3942, -1.1950, -0.1545,
  //	     -0.1360, 0.4526, -0.1936);
  TF2 *f2(0), *f2b(0);
  TF1 *f1(0); 
  double norm_mip = 0.045; // 0.065;
  double ptmin = (pt>0 ? min(pt,20.) : 3.);//8.);//20.);
  int ipt0 = (pt>0 ? h3->GetXaxis()->FindBin(pt) : 0);
  int iptmin = (pt>0 ? ipt0 : h3->GetXaxis()->FindBin(ptmin));
  //double ipt1 = h3->GetXaxis()->FindBin(49.);
  double ipt1 = h3->GetXaxis()->FindBin(199.);
  bool goodfit = false;
  for (int ipt = h3->GetNbinsZ(); ipt != iptmin-1; --ipt) {

    if (ipt0!=0 && ipt!=ipt0) continue;
    double pt1 = h3->GetXaxis()->GetBinLowEdge(ipt);
    double pt2 = h3->GetXaxis()->GetBinLowEdge(ipt+1);
    double pt = 0.5*(pt1+pt2);
    double dpt = 0.5*(pt2-pt1);

    //if (pt2>50) continue; // VLE sample
    if (pt2>200) continue; // VLE+ME sample

    TH2D *h2 = (TH2D*)h3->Project3D("yz");
    h2->GetYaxis()->SetRangeUser(0,1.5);

    TProfile *px = h2->ProfileX();
    TProfile *py = h2->ProfileY();
    //px->Draw("SAMEP");
    //TGraphErrors *gy = new TGraphErrors(py->ProjectionX());
    //for (int i = 0; i != gy->GetN(); ++i) {
    //gy->SetPoint(i, py->GetBinContent(i), py->GetBinCenter(i));
    //}
    //gy->Draw("SAMEP");

    //double minH = minpt_hcal/pt1;
    //double minE = mip_ecal/pt1;
    double minH0 = minpt_hcal/pt1;
    int iminH0 = px->FindBin(minH0);
    double minH = px->GetBinLowEdge(iminH0+1);
    //double minE = py->GetBinLowEdge(py->FindBin(minpt_ecal/pt1)+1);
    double minE0 = mip_ecal/pt1;
    int iminE0 = py->FindBin(minE0);
    double minE = py->GetBinLowEdge(iminE0+1);

    // Select just one pT bin
    for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	
	// normalize E=0 bin, which is due to MIPS
	double ET = h2->GetYaxis()->GetBinLowEdge(j);
	//double norm = (ET<mip_ecal/pt1 ? 0.5 : 1);
	//double norm = (ET<mip_ecal/pt1 ? 0.065 : 1);
	//double norm = (ET<mip_ecal/pt1 ? norm_mip : 1);
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
    //double r0 = (h2->GetMean(1)+h2->GetMean(2));
    // double s0 = sqrt(pow(h2->GetRMS(1),2) + pow(h2->GetRMS(2),2));
    double r0 = (hE->GetMean()+hH->GetMean());
    double s0 = sqrt(pow(hE->GetRMS(),2) + pow(hH->GetRMS(),2));

    TH1D *h = new TH1D("h",";H / T (GeV);E / T (GeV);",75,0,1.5);//100,0,2);
    h->SetMinimum(0);
    h->SetMaximum(1.5);//2);
    
    extraText = "Simulation";
    lumi_8TeV = "Pythia 6 Z2*";
    TCanvas *c1 = tdrCanvas("c1", h,2,11, true);
  
    gStyle->SetPalette(1);
    tdrDraw(h2,"COLZ");
    //tdrDraw(h2n,"COLZ");

    TBox *box = new TBox();
    box->SetFillStyle(3001);
    box->SetFillColor(kGray);
    //box->DrawBox(0.00,4./pt2,1.25,15./pt1);
    

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045);
    tex->DrawLatex(0.45,0.87,Form("%1.1f < T < %1.1f GeV",pt1,pt2));
    //tex->DrawLatex(0.45,0.82,Form("Pion gun, mag. field off",pt1,pt2));
    tex->DrawLatex(0.45,0.82,Form("Pion gun, mag. field on",pt1,pt2));
    tex->DrawLatex(0.22,0.65,"E-hadrons");
    tex->DrawLatex(0.46,0.44,"EH-hadrons");
    tex->DrawLatex(0.71,0.22,"H-hadrons");
    tex->SetTextColor(kGray+2); tex->SetTextSize(0.035);
    //tex->DrawLatex(0.18,0.15,"Lost 4<E<15 (bug)");
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.75,0.19,Form("(#times1/%d)",int(1./norm_mip)));

    TArrow *arr = new TArrow(); arr->SetNDC(kTRUE);
    arr->DrawArrow(0.1,1.0,0.05,1.0,0.03);//0.21,0.72,0.18,0.72);
    arr->DrawArrow(0.57,0.57,0.52,0.52,0.03);//0.47,0.47,0.43,0.43);
    arr->DrawArrow(1.04,0.13,1.01,0.05,0.03);//0.75,0.24,0.70,0.18);

    TLine *l = new TLine();
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1,1,0);
    l->DrawLine(0,0.667,0.667,0);
    l->DrawLine(0,0.333,0.333,0);
    l->SetLineStyle(kSolid);
    l->DrawLine(minH,1-minH,1-minE,minE);
    l->DrawLine(minH,minE,minH,1.2);//1-minH);
    l->DrawLine(minH,minE,1.4,minE);//1-minE,minE);
    l->SetLineStyle(kDashed);
    l->DrawLine(minH0,minE0,minH0,1-minH0);
    l->DrawLine(minH0,minE0,1-minE0,minE0);

    double ptx = 45;
    // Before photon fix
    //double mine = (pt1<ptx ? 1./pt1 : 15.5/pt1);
    //double minh = (pt1<ptx ? 0.8/pt1 : 0.1);
    //double maxe = (pt1<ptx ? 1.0-minh : 0.8);
    //double maxh = (pt1<ptx ? 1.0-mine : 0.8);
    double mine = max(0.0,minE);//mip_ecal/pt1;
    double minh = max(0.0,minH);//max(0.1,minpt_hcal/pt1);
    double maxe = min(1.0,1.0-minh);
    double maxh = min(1.0,1.0-mine);
    //double maxe = min(0.65,1.0-minh);
    //double maxh = min(0.80,1.0-mine);
    // Instantiate these only once so we can recycle fit from previous round
    f2 = (f2 ? f2 : new TF2("f2",func,minh,maxh,mine,maxe,6));
    f2b = (f2b ? f2b : new TF2("f2b",func,0.,1,0.,1.,6));
    //f2->SetParameters(6.1, 0.858, 0.0738, 0.812, 0.402, 0.356);
    //f2->SetParameters(6.1, 0.70, 0.12, 0.58, 0.1, 0.5);
    if (!goodfit) 
      f2->SetParameters(3.14, 0.866, 0.0970, 0.638, 0.577, 0.351); // 48-50 GeV
    if (pt1<25) {
      f2->FixParameter(4, 1.); // low pT showers all in ECAL
      //f2->FixParameter(5, 1.); // keep some constraint on E/H scatter
    }
    
    //feh->SetParameter(0,1); feh->SetParameter(1,0);
    //double a0 = feh->Eval(pt1);
    //feh->SetParameter(0,0); feh->SetParameter(1,1);
    //double b0 = feh->Eval(pt1);
    //f2->SetParameters(6.1, a0, 0.12, b0, 0.1, 0.5);

    //f2->SetParLimits(4,0.35,0.75);
    //if (pt1<ptx) f2->FixParameter(4,0.5);
    ///if (pt1<ptx) f2->FixParameter(5,1.0);
    h2->Fit(f2,ipt0 ? "RN" : "QRN");
    h2->Fit(f2,ipt0 ? "RNME" : "QRNME");
    //f2->Draw("contl same");
    f2->SetContour(4);
    f2b->SetContour(4);
    for (int i=0; i!=f2->GetNpar(); ++i)
      f2b->SetParameter(i,f2->GetParameter(i));
    f2b->SetLineColor(kBlue+1);
    f2b->Draw("SAME");
    f2->Draw("SAME");
    
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
    l->SetLineWidth(2);
    l->SetLineStyle(kSolid);
    l->SetLineColor(kBlue+2);
    l->DrawLine(0,re,rh,0);
    TEllipse *tel = new TEllipse();
    tel->SetFillStyle(kNone);
    //tel->DrawEllipse(hf,0,0.02,0.02,0,360,0);
    //tel->DrawEllipse(0,(1-hf),0.02,0.02,0,360,0);
    tel->DrawEllipse(1-ef,ef,0.02,0.02,0,360,0);
    tel->SetLineColor(kMagenta);
    tel->DrawEllipse(rh*(1-ef),re*ef,0.02,0.02,0,360,0);
				      
    gPad->RedrawAxis();
    
    if (ipt0) {
      c1->SaveAs(Form("pdf/pfhadrons_2D_%1.1f-%1.1f.pdf",pt1,pt2));
    }
    else {

      // Animated gif scanning from high pT to low pT
      // Remove gif before starting to make sure it won't get messed up?
      if (pt1>10 && ipt<=ipt1)
	c1->SaveAs(Form("pdf/pfhadrons_2D_anim.gif+%d",(ipt1-ipt)*100));

      //{
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

      delete c1;
      delete h;
      delete h2;
      delete px;
      delete py;
      //delete f2;
      //delete f2b;
    }
  } // for ipt

  if (ipt0==0) {

    TH1D *h = new TH1D("h",";T (GeV);Hadron response;",200,0,200);
    h->GetXaxis()->SetRangeUser(ptmin,200);
    h->SetMinimum(0);
    h->SetMaximum(1.1);
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();

    TCanvas *c1 = tdrCanvas("c1",h,2,11,true);
    gPad->SetLogx();

    tdrDraw(ghr, "P", kFullCircle, kRed);
    //tdrDraw(gar, "P", kFullDiamond, kBlack);
    tdrDraw(ger, "P", kFullSquare, kBlue);
    //tdrDraw(gh4, "P", kOpenCircle, kRed);

    tdrDraw(grh, "P", kOpenTriangleDown, kMagenta+1);//kRed);
    tdrDraw(gre, "P", kOpenTriangleUp, kCyan+1);//kBlue);
    grr->SetMarkerSize(2.0);
    //tdrDraw(grr, "P", kFullStar, kMagenta+1);
    grd->SetMarkerSize(2.0);
    //tdrDraw(grd, "P", kOpenStar, kMagenta+1);

    //TF1 *fhr = new TF1("fhr","[2]-[0]*pow(x,-[1])",ptminfit,ptmaxfit);
    //fhr->SetParameters(1,0.15,1);
    // Use same function as in 1D fit
    TF1 *fhr = new TF1("fhr","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",4,200);
    fhr->SetParameters(1.0806,-0.5,-0.15,-0.3);
    fhr->FixParameter(0, 1.0806); // constrain to old VLE+ME sample
    fhr->SetParLimits(2,-0.18,-0.13);
    cout << "*** Fitting fhr ***" << endl;
    ghr->Fit(fhr,"RN");
    fhr->DrawClone("SAME");
    fhr->SetLineStyle(kDashed);
    fhr->SetRange(-fhr->GetParameter(3),200);
    fhr->DrawClone("SAME");

    TF1 *frh = new TF1("frh","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    frh->SetParameters(0.9*fhr->GetParameter(0), fhr->GetParameter(1),
		       fhr->GetParameter(2), fhr->GetParameter(3));
    frh->FixParameter(1,fhr->GetParameter(1));
    frh->SetParLimits(2,-0.18,-0.13);
    frh->SetLineColor(kMagenta+1);//Red);
    cout << "*** Fitting frh ***" << endl;
    grh->Fit(frh,"RN");
    frh->DrawClone("SAME");
    frh->SetLineStyle(kDashed);
    frh->SetRange(-frh->GetParameter(3),200);
    frh->DrawClone("SAME");

    TF1 *fr = new TF1("fr","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    fr->SetParameters(0.95*frh->GetParameter(0), frh->GetParameter(1),
		      frh->GetParameter(2), frh->GetParameter(3));
    fr->SetParLimits(2,-0.18,-0.13);
    fr->SetLineColor(kMagenta+1);
    cout << "*** Fitting fr ***" << endl;
    grr->Fit(fr,"RN");
    //fr->DrawClone("SAME");
    fr->SetLineStyle(kDashed);
    fr->SetRange(-fr->GetParameter(3),200);
    //fr->DrawClone("SAME");

    TF1 *fre = new TF1("fre","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
    fre->SetParameters(0.95*fr->GetParameter(0), fr->GetParameter(1),
		       fr->GetParameter(2), fr->GetParameter(3));
    fre->SetParLimits(2,-0.18,-0.13);
    fre->SetLineColor(kCyan+1);
    cout << "*** Fitting fre ***" << endl;
    gre->Fit(fre,"RN");
    fre->DrawClone("SAME");
    fre->SetLineStyle(kDashed);
    fre->SetRange(-fre->GetParameter(3),200);
    fre->DrawClone("SAME");

    TF1 *fer = new TF1("fer","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",3,200);
    fer->SetParameters(0.98*fre->GetParameter(0), fre->GetParameter(1),
		       fre->GetParameter(2), fre->GetParameter(3));
    fre->SetParLimits(2,-0.18,-0.13);
    fer->SetLineColor(kBlue);
    cout << "*** Fitting fer ***" << endl;
    ger->Fit(fer,"RN");
    fer->DrawClone("SAME");
    fer->SetLineStyle(kDashed);
    fer->SetRange(-fer->GetParameter(3),200);
    fer->DrawClone("SAME");

    TF1 *fd = new TF1("fd","[0]+[1]*pow(x,[2])",9,200);
    fd->SetParameters(0.01,-0.01,-0.15);
    fd->SetLineColor(kMagenta+1);
    cout << "*** Fitting fd ***" << endl;
    grd->Fit(fd,"RN");
    //fd->DrawClone("SAME");
    fd->SetLineStyle(kDashed);
    fd->SetRange(3,200);
    //fd->DrawClone("SAME");

    //TLegend *leg = tdrLeg(0.55,0.22,0.85,0.52);
    TLegend *leg = tdrLeg(0.55,0.22,0.85,0.42);
    leg->AddEntry(ghr,"H-hadron","P");
    leg->AddEntry(grh,"EH-had. (HCAL)","P");
    //leg->AddEntry(grr,"0.5#times(R_{E}+R_{H})","P");
    leg->AddEntry(gre,"EH-had. (ECAL)","P");
    leg->AddEntry(ger,"E-hadron","P");
    //leg->AddEntry(grd,"0.5#times(R_{H}-R_{E})","P");


    TH1D *h2 = (TH1D*)h->Clone("h2");
    h2->GetYaxis()->SetTitle("Hadron resolution (#sigma / #mu)");

    TCanvas *c2 = tdrCanvas("c2",h2,2,33,true);
    gPad->SetLogx();

    tdrDraw(ghs, "P", kFullCircle, kRed);
    //tdrDraw(gas, "P", kFullDiamond, kBlack);
    tdrDraw(ges, "P", kFullSquare, kBlue);
    grs->SetMarkerSize(2.0);
    tdrDraw(grs, "P", kFullStar, kMagenta+1);

    TF1 *fhs = new TF1("fhs","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",7,200);
    fhs->SetParameters(1,1,0.1);
    //fhs->FixParameter(0,0);
    fhs->SetLineColor(kRed);
    cout << "*** Fitting fhs ***" << endl;
    ghs->Fit(fhs,"RN");
    fhs->DrawClone("SAME");
    fhs->SetLineStyle(kDashed);
    fhs->SetRange(3,200);
    fhs->DrawClone("SAME");

    TF1 *fes = new TF1("fes","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",7,200);
    fes->SetParameters(1,1,0.1);
    //fes->FixParameter(0,0);
    fes->SetLineColor(kBlue);
    cout << "*** Fitting fes ***" << endl;
    ges->Fit(fes,"RN");
    fes->DrawClone("SAME");
    fes->SetLineStyle(kDashed);
    fes->SetRange(3,200);
    fes->DrawClone("SAME");

    TF1 *fs = new TF1("fs","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",7,200);
    fs->SetParameters(1,1,0.1);
    //fs->FixParameter(0,0);
    fs->SetLineColor(kMagenta+1);
    cout << "*** Fitting fs ***" << endl;
    grs->Fit(fs,"RN");
    fs->DrawClone("SAME");
    fs->SetLineStyle(kDashed);
    fs->SetRange(3,200);
    fs->DrawClone("SAME");
    
    TLegend *leg2 = tdrLeg(0.55,0.45,0.85,0.60);
    leg2->AddEntry(ges,"E-hadron","P");
    leg2->AddEntry(grs,"EH-hadron","P");
    leg2->AddEntry(ghs,"H-hadron","P");


    TH1D *h3 = (TH1D*)h->Clone("h3");
    h3->GetYaxis()->SetTitle("EH-hadron ECAL deposition");

    TCanvas *c3 = tdrCanvas("c3",h3,2,11,true);
    gPad->SetLogx();

    gfm->SetMarkerStyle(kOpenDiamond);
    gfm->Draw("SAMEP");

    TLegend *leg3 = tdrLeg(0.55,0.30,0.85,0.60);
    leg3->AddEntry(gfm,"Centroid","P");

    cout << endl;
    cout << "New feh: r+/-dr" << endl;
    cout <<Form("feh->SetParameters(0.5, 0.5,\n"
		"                   %1.4f, %1.4f, %1.4f, %1.2f,\n"
		"                   %1.4f, %1.4f, %1.4f, %1.2f);\n",
		frh->GetParameter(0), frh->GetParameter(1),
		frh->GetParameter(2), frh->GetParameter(3),
		fre->GetParameter(0), fre->GetParameter(1),
		fre->GetParameter(2), fre->GetParameter(2));
	       //fr->GetParameter(0), fr->GetParameter(1), fr->GetParameter(2),
	       //fd->GetParameter(0), fd->GetParameter(1), fd->GetParameter(2));
    cout << endl;

    c1->SaveAs("pdf/pfhadrons_2Dfit_R.pdf");
    c2->SaveAs("pdf/pfhadrons_2Dfit_S.pdf");
    c3->SaveAs("pdf/pfhadrons_2Dfit_F.pdf");
  }

} // draw2D


void pfhadrons::drawJER() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile("pfhadrons_jer.root","READ");
  assert(f && !f->IsZombie());

  vector<string> modes;
  modes.push_back("Ec");
  modes.push_back("Hec");
  modes.push_back("EHc");
  modes.push_back("Ac");
  
  int color[] = {kBlue, kRed, kMagenta+1, kBlack};
  int marker[] = {kFullSquare, kFullCircle, kFullStar, kOpenCircle};

  vector<TGraphErrors*> vs(modes.size());
  vector<TGraphErrors*> vr(modes.size());
  for (unsigned int i = 0; i != modes.size(); ++i) {
    vs[i] = (TGraphErrors*)f->Get(Form("jerc_%s",modes[i].c_str()));
    assert(vs[i]);
    vr[i] = (TGraphErrors*)f->Get(Form("jesc_%s",modes[i].c_str()));
    assert(vr[i]);
  }
  
  TH1D *h1 = new TH1D("h1",";T (GeV);Hadron resolution",1990,1,200);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->SetMaximum(1.1);

  TCanvas *c1 = tdrCanvas("c1",h1 ,2,33,true);
  c1->SetLogx();

  for (unsigned int i = 0; i != vs.size(); ++i) {
    tdrDraw(vs[i],"P",marker[i],color[i]);
  }
  c1->SaveAs("pdf/pfhadrons_jer.pdf");
  
  TH1D *h2 = new TH1D("h2",";T (GeV);Residual hadron scale",1990,1,200);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetMaximum(1.6);
  h2->SetMinimum(0.7);
  
  TCanvas *c2 = tdrCanvas("c2",h2,2,11,true);
  c2->SetLogx();
  
  for (unsigned int i = 0; i != vr.size(); ++i) {
    tdrDraw(vr[i],"P",marker[i],color[i]);
  }
  c2->SaveAs("pdf/pfhadrons_jes.pdf");
  
}
