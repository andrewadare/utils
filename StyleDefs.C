// To use one of these styles:
// gROOT->Macro("StyleDefs.C");
// gROOT->SetStyle("custom1")
//
// To reset style to the default:
// gROOT->SetStyle("Modern")
//
// Palettes: 
// use gStyle->SetPalette(ncolors,0) with ncolors = 
// 51 blue, 52 grey, 53 "radiator", 54 blue/gray/yellow
// 55 (rainbow) is default.

#include "TStyle.h"

TStyle* s1 = 0;
TStyle* s2 = 0;
int niceBlue[100];
int niceRed[100];

void StyleDefs() {

  {
    // niceBlue: 4 colors, 100 steps  
    Double_t r[]  = { 1.00, 0.50, 0.00, 0.05};
    Double_t g[]  = { 1.00, 0.90, 0.20, 0.15};
    Double_t b[]  = { 1.00, 1.00, 1.00, 0.50};
    Double_t upto[] = { 0.0, 0.10, 0.80, 1.00 };
    int ct = TColor::CreateGradientColorTable(4,upto,r,g,b,100);
    for (int i=0; i<100; i++) niceBlue[i] = ct+i;
  }
  {
    // niceRed: 4 colors, 100 steps  
    Double_t r[]  = { 1.00, 1.00, 1.00, 0.50};
    Double_t g[]  = { 1.00, 0.90, 0.20, 0.15};
    Double_t b[]  = { 1.00, 0.50, 0.00, 0.10};
    Double_t upto[] = { 0.0, 0.10, 0.50, 1.00 };
    int ct = TColor::CreateGradientColorTable(4,upto,r,g,b,100);
    for (int i=0; i<100; i++) niceRed[i] = ct+i;
  }


  // "s1" style ===================================================
  s1 = new TStyle("custom1", "custom style 1");

  s1->SetPalette(100,niceBlue);

  // TLatex
  s1->SetTextFont(42);

  // Canvas and pad
  s1->SetCanvasBorderMode(0);
  s1->SetCanvasColor(0);
  s1->SetPadBorderMode(0);
  s1->SetPadColor(0);
  
  // Margins
  s1->SetPadBottomMargin(0.14);
  s1->SetPadLeftMargin(0.14);
  s1->SetPadRightMargin(0.05);
  
  // Frame
  s1->SetFrameBorderMode(0);
  s1->SetFrameFillColor(0);
  
  // Axes
  s1->SetPadTickX(0);
  s1->SetPadTickY(1);
  s1->SetTitleOffset(1.3,"x");
  s1->SetTitleOffset(1.4,"y");
  s1->SetTitleOffset(1.4,"z");
  
  // Axis label/title sizes
  s1->SetLabelSize(0.045,"x");
  s1->SetTitleSize(0.045,"x");
  s1->SetLabelSize(0.045,"y");
  s1->SetTitleSize(0.045,"y");
  s1->SetLabelSize(0.045,"z");
  s1->SetTitleSize(0.045,"z");

  // Title/axis fonts
  s1->SetLabelFont(42,"x");
  s1->SetTitleFont(42,"x");
  s1->SetLabelFont(42,"y");
  s1->SetTitleFont(42,"y");
  s1->SetLabelFont(42,"z");
  s1->SetTitleFont(42,"z");

  // Title
  s1->SetTitleSize(0.050,"");
  s1->SetTitleFont(42,"");
  s1->SetTitleAlign(23);
  s1->SetTitleX(0.5);
  s1->SetTitleBorderSize(0);
  s1->SetTitleFillColor(0);
  s1->SetTitleStyle(0);

  // Stat box turned off
  s1->SetOptStat(0);
  s1->SetStatY(0.935);
  s1->SetStatColor(0);
  s1->SetStatBorderSize(1);
  s1->SetStatFont(42);

  // Histogram options
  s1->SetHistLineColor(kBlue+2);
  s1->SetHistLineWidth(2);
  s1->SetHistTopMargin(0.10); // Default 0.05

  // TF1 options
  s1->SetFuncWidth(2);
  s1->SetFuncColor(kRed);

  // Legend options
  s1->SetLegendBorderSize(1);
  s1->SetLegendFillColor(0);
  s1->SetLegendFont(42);

  // Add new style defs below ========================================
  //
  // "s2" style ===================================================
  s2 = new TStyle("custom2", "custom style 2");

  s2->SetPalette(51,0);

  // Canvas and pad
  s2->SetCanvasBorderMode(0);
  s2->SetCanvasColor(0);
  s2->SetPadBorderMode(0);
  s2->SetPadColor(0);
  
  // Margins
  s2->SetPadBottomMargin(0.14);
  s2->SetPadLeftMargin(0.14);
  s2->SetPadRightMargin(0.14);
  
  // Frame
  s2->SetFrameBorderMode(0);
  s2->SetFrameFillColor(0);
  
  // Axes
  s2->SetPadTickX(0);
  s2->SetPadTickY(1);
  s2->SetTitleOffset(1.3,"x");
  s2->SetTitleOffset(1.4,"y");
  s2->SetTitleOffset(1.4,"z");
  
  // Axis label/title sizes
  s2->SetLabelSize(0.045,"x");
  s2->SetTitleSize(0.045,"x");
  s2->SetLabelSize(0.045,"y");
  s2->SetTitleSize(0.045,"y");
  s2->SetLabelSize(0.045,"z");
  s2->SetTitleSize(0.045,"z");

  // Title/axis fonts
  s2->SetLabelFont(42,"x");
  s2->SetTitleFont(42,"x");
  s2->SetLabelFont(42,"y");
  s2->SetTitleFont(42,"y");
  s2->SetLabelFont(42,"z");
  s2->SetTitleFont(42,"z");

  // Title
  s2->SetTitleSize(0.050,"");
  s2->SetTitleFont(42,"");
  s2->SetTitleAlign(23);
  s2->SetTitleX(0.5);
  s2->SetTitleBorderSize(0);
  s2->SetTitleFillColor(0);
  s2->SetTitleStyle(0);

  // Stat box turned off
  s2->SetOptStat(0);
  s2->SetStatY(0.935);
  s2->SetStatColor(0);
  s2->SetStatBorderSize(1);
  s2->SetStatFont(42);

  // Histogram options
  s2->SetHistLineColor(kBlue+2);
  s2->SetHistLineWidth(2);
  s2->SetHistTopMargin(0.10); // Default 0.05

  // TF1 options
  s2->SetFuncWidth(2);
  s2->SetFuncColor(kRed);

  // Legend options
  s2->SetLegendBorderSize(1);
  s2->SetLegendFillColor(0);
  s2->SetLegendFont(42);
}
