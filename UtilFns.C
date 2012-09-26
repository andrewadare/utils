#include "TCanvas.h"
#include "TFile.h"
#include "TList.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TString.h"
#include "TLatex.h"
#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TSystem.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphTime.h"

// --- Functions ---
TCanvas* DrawObject(TObject* obj, 
		    TString drawopt="", 
		    TString title="",
		    TObjArray* cList=0,
		    double xpx=700, double ypx=500);

void SaveCanvases(TObjArray* canvases, const char* fileName);
TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir = "");
void SaveCanvasesFromFile(const char* rootFile, 
			  const char* targetDir, 
			  const char* tag, 
			  const char* fileType);
int PrintPDF(TObjArray* cList, TString base, TString opt = "pdf");
void SetHistProps(TH1* h,
		  Int_t linecolor,
		  Int_t fillcolor,
		  Int_t markercolor,
		  Int_t markerstyle,
		  Double_t markersize);

TGraphTime* Animation(TObjArray* moveObjs, TObjArray* statObjs,
		      TString opt="", int sleeptime=50);

void SaveCanvases(TObjArray* canvases, const char* fileName)
{
  TFile* f = new TFile(fileName, "recreate");

  if (!canvases)
    gROOT->Error("UtilFns::SaveCanvases()", "!canvases");

  for (int n=0; n<canvases->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)canvases->At(n);
    if (c) {
      c->Write(c->GetTitle());
    }
    else
      gROOT->Warning("UtilFns::SaveCanvases()", "!c %d", n);
  }
  if (1)
    gROOT->Info("", "Wrote %s", f->GetName());

  f->Close();
  return;
}

void SaveCanvasesFromFile(const char* rootFile, 
			  const char* targetDir, 
			  const char* tag, 
			  const char* fileType)
{
  // Get a list of canvases from rootFile into array, then save each
  // to its own file in targetDir/. fileType = "eps", "pdf", "C",
  // "png", etc. Not all formats have been tested.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString name = "";
  TString base(targetDir);
  TFile *cFile = new TFile(rootFile, "read");
  TObjArray* cList = GetObjectsFromFile(*cFile, "TCanvas");

  if (!cList) {
    gROOT->Error("UtilFns::SaveCanvasesFromFile()", "!cList");
  }
  
  for (int n=0; n<cList->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)cList->At(n);
    if (c) {
      name = "";
      name = base;
      name += TString("/");
      name += TString(fileType);
      name += TString("/");
      name += TString(c->GetTitle());
      name += TString(".");
      name += TString(fileType);
      
      if (0)
	gROOT->Info("", "%s", name.Data());

      c->Draw();
      c->Modified();
      c->Update();
      c->SaveAs(name.Data());
    }
    else
      gROOT->Error("SaveCanvasesFromFile()", "!c");
  }

  if (1) {
    PrintPDF(cList, Form("%s/all-figs%s", targetDir, tag), "pdf");
  }
  
  return;
}

TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir)
{
  file.cd(dir.Data());

  TObjArray* objList = new TObjArray();
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  
  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    if (0) 
      printf("%10s %20s\n", className.Data(), keyName.Data());
    
    if (className.Contains(clname)) {
      objList->Add(gDirectory->Get(keyName.Data()));
    }
  }

  return objList;
}

int PrintPDFs(TObjArray* cList, TString dir="./", TString opt="")
{
  Int_t nPrinted = 0;  
  TString ext = ".pdf";

  if (opt.Contains("ps"))
    ext = ".ps";
  if (opt.Contains("eps"))
    ext = ".eps";
  
  if (!cList) {
    gROOT->Error("PrintPDF()", "no cList!");
    return -1;
  }

  TCanvas* c = 0;
  for (int i=0; i<cList->GetEntries(); i++) {
    TObject* obj = cList->At(i);
    if (TString(obj->ClassName()).Contains("TCanvas"))
      c = (TCanvas*)obj;
    else {
      gROOT->Warning("PrintPDFs()", 
		     "list contains non-canvas object");
      continue;
    }
    if (!c) {
      gROOT->Error("PrintPDFs()", "no canvas!");
      return -1;
    }

    TString fileName = TString(c->GetName()) + ext;
    if (gSystem->FindFile(dir.Data(), fileName))
      gROOT->Info("PrintPDFs()", "Overwriting %s", fileName.Data());
    else    
      fileName = dir + fileName;

    c->Print(fileName.Data(), ext);
    nPrinted++;
  }
  return nPrinted;
}

int PrintPDF(TObjArray* cList, TString base, TString opt)
{
  TLatex ltx;
  ltx.SetNDC();
  Int_t nPrinted = 0;  
  TString psOut("");
  
  TString ext = ".pdf";

  if (opt.Contains("ps"))
    ext = ".ps";
  if (opt.Contains("eps"))
    ext = ".eps";

  if (!cList) {
    gROOT->Error("PrintPDF()", "no cList!");
    return -1;
  }

  TCanvas* c = 0;
  for (int i=0; i<cList->GetEntries(); i++) {
    TObject* obj = cList->At(i);
    if (TString(obj->ClassName()).Contains("TCanvas"))
      c = (TCanvas*)obj; //(TCanvas*)cList->At(i); 
    else {
      gROOT->Warning("PlotUtils::print_pdf()", "list contains non-canvas object");
      continue;
    }
    if (!c) {
      gROOT->Error("PlotUtils::print_pdf()", "no canvas!");
      return -1;
    }

    // Slide numbering
    if (opt.Contains("number")) {
      c->cd();
      ltx.DrawLatex(0.95, 0.95, Form("%d", nPrinted+1));
    }
    
    // Multipage ps and pdf files
    if (nPrinted==0) {
      psOut = base + ext + "[";
      c->Print(psOut.Data()); // opens ps but doesn't print it
    }
    if (i < cList->GetEntries()) {
      psOut = base + ext;
      if(ext.Contains("pdf")) {
	TString pageName = Form("Title:%s", c->GetTitle());
	c->Print(psOut.Data(), pageName.Data());
      }
      else
	c->Print(psOut.Data());
    }
    if (i==cList->GetEntries()-1) {
      psOut = base + ext + "]";
      c->Print(psOut.Data()); // closes ps but doesn't print it
    }

    // Also print pages as individual files, if requested
    if (opt.Contains("singles")) {
      TString dir = gSystem->DirName(base.Data());
      TString fileName = dir + "/" + TString(c->GetName()) + ext;
      c->Print(fileName.Data());
    }

    nPrinted++;
  }
    
  if (ext.Contains("ps")) {
    TString cmd = "ps2pdf " + base + ext + " " + base + ".pdf";
    gSystem->Exec(cmd.Data()); 
  }

  return 0;
}

TCanvas* DrawObject(TObject* obj, 
		    TString drawopt, 
		    TString title,
		    TObjArray* cList,
		    double xpx, double ypx)
{
  // Draw a TH1, TGraph, or anything with a Draw() method, in a new canvas.
  // Use opt to set drawing options.

  static int ci = 0;
  double x = xpx > 0 ? xpx : 700;
  double y = ypx > 0 ? ypx : 500;
  TCanvas* c = new TCanvas(Form("c%d",ci),Form("c%d",ci), x, y);
  ci++;

  if (!title.IsNull() && obj->InheritsFrom("TNamed")) {
    dynamic_cast<TNamed*>(obj)->SetTitle(title.Data());
    c->SetTitle(title.Data());
  }

  if (drawopt.Contains("clone")) {
    drawopt.ReplaceAll("clone", ""); // clear unwanted (c,l,e) options
    obj->DrawClone(drawopt.Data());
  }
  else
    obj->Draw(drawopt.Data());
  
  if (cList) 
    cList->Add(c);

  return c;
}

void SetHistProps(TH1* h,
		  Int_t linecolor,
		  Int_t fillcolor,
		  Int_t markercolor,
		  Int_t markerstyle,
		  Double_t markersize) 
{
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
}

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize) 
{
  g->SetLineColor(linecolor);
  g->SetMarkerColor(markercolor);
  g->SetMarkerStyle(markerstyle);
  g->SetMarkerSize(markersize);
  g->SetLineWidth(2);
}

TGraphTime* Animation(TObjArray* moveObjs, TObjArray* statObjs,
		      TString opt, int sleeptime)
{
  static int iAnim = 0; iAnim++;
  int nFrames = moveObjs->GetEntries();
  gROOT->Info("Animation()", "Creating %d frame sequence...", nFrames);
  TGraphTime* anim = new TGraphTime(nFrames,0,0,1,1);
  anim->SetName(Form("anim%d", iAnim));
  
  for (int n=0; n<nFrames; n++) {
    
    // Add stationary objects to this frame
    if (statObjs) {
      for (int i=0; i<statObjs->GetEntries(); i++) {
	TObject* sta = statObjs->At(i);
	TString drawOpt = i ? "same" : "";
	if (!opt.IsNull())
	  drawOpt += opt;
	anim->Add(sta, n, drawOpt);
      }
    }

    // Add changing objects
    TObject* mov = moveObjs->At(n);
    TString drawOpt2 = statObjs ? "same" : "";
    if (!opt.IsNull())
      drawOpt2 += opt;
    anim->Add(mov, n, drawOpt2);

  }

  anim->SetTitle("animation");
  anim->SetSleepTime(sleeptime); // ms (default = 0)
  return anim;
}
