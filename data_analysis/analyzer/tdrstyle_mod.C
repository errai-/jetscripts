#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"

// tdrGrid: Turns the grid lines on (true) or off (false)

void tdrGrid(bool gridOn) {
  TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); assert(tdrStyle);
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis
void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.15);//0.13);
  tdrStyle->SetPadLeftMargin(0.15);//0.16);
  tdrStyle->SetPadRightMargin(0.05);//0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.1);//0.9);
  tdrStyle->SetTitleYOffset(1.25); // => 1.15 if exponents
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  // Additional settings for QCD-10-011
  tdrStyle->SetLegendBorderSize(0);

  tdrStyle->cd();

}


void cmsPrel4x4(double intLumi=-1, bool wide = true) {

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  
  latex->SetTextAlign(31); // align right
  //latex->DrawLatex(wide ? 0.98 : 0.95, 0.96, "#sqrt{s} = 7 TeV");
  latex->DrawLatex(wide ? 0.97 : 0.95, 0.96, "#sqrt{s} = 8 TeV");
  if (intLumi > 0.) {
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(wide ? 0.08 : 0.15, 0.96,
		     //Form("CMS preliminary, L = %.2g pb^{-1}",intLumi));
		     //Form("CMS preliminary, L = %.0f pb^{-1}",intLumi));
		     Form("CMS preliminary, L = %.2g fb^{-1}",intLumi*0.001));
  }
  else if (intLumi==0) { // simulation
    latex->SetTextAlign(11); // align left
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Fall10 QCD)");
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Summer11 QCD)");
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Summer11)");
    latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation 2014");
  }
  else {
    latex->SetTextAlign(11); // align left
    //latex->DrawLatex(0.15,0.96,"CMS preliminary 2010");
    //latex->DrawLatex(0.15,0.96,"CMS preliminary 2011");
    latex->DrawLatex(0.15,0.96,"CMS preliminary 2014");
  }
} // cmsPrel


void cmsPrel(double intLumi=-1, bool wide = false) {

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  
  latex->SetTextAlign(31); // align right
  //latex->DrawLatex(wide ? 0.98 : 0.95, 0.96, "#sqrt{s} = 7 TeV");
  latex->DrawLatex(wide ? 0.98 : 0.95, 0.96, "#sqrt{s} = 8 TeV");
  if (intLumi > 0.) {
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(wide ? 0.06 : 0.15, 0.96,
		     //Form("CMS preliminary, L = %.2g pb^{-1}",intLumi));
		     //Form("CMS preliminary, L = %.0f pb^{-1}",intLumi));
		     Form("CMS preliminary, L = %.2g fb^{-1}",intLumi*0.001));
  }
  else if (intLumi==0) { // simulation
    latex->SetTextAlign(11); // align left
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Fall10 QCD)");
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Summer11 QCD)");
    //latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Summer11)");
    latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation 2014");
  }
  else {
    latex->SetTextAlign(11); // align left
    //latex->DrawLatex(0.15,0.96,"CMS preliminary 2010");
    //latex->DrawLatex(0.15,0.96,"CMS preliminary 2011");
    latex->DrawLatex(0.15,0.96,"CMS preliminary 2014");
  }
} // cmsPrel

void cmsFinal(double intLumi=-1, bool wide = false) {

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(wide ? 0.98 : 0.95, 0.96, "#sqrt{s} = 7 TeV");
  if (intLumi > 0.) {
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(wide ? 0.06 : 0.15, 0.96,
		     //Form("CMS, L = %.2g pb^{-1}",intLumi)); // JINST
		     Form("CMS, L = %.2g fb^{-1}",intLumi*0.001)); // JINST
  }
  else if (intLumi==0) { // simulation
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(wide ? 0.06 : 0.15, 0.96, "CMS simulation (Pythia Z2)");
  }
  else {
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.15,0.96,"CMS 2011");
  }
} // cmsPrel

//cmsPrel(); // to print just CMS and \sqrt{s}
//cmsPrel(400);  // to print also the integrated luminosity.
