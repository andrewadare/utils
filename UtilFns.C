// UtilFns.C
//
// ROOT functions for plotting, image output, animation, etc.
// Get the latest version:
// git clone https://github.com/andrewadare/utils.git
//
// To use in your interpreted macro, add this line:
// gROOT->LoadMacro("/path/to/UtilFns.C");
//
// Usage in compiled code may require some hacking.
// -Andrew Adare 4/24/2013

#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphTime.h"
#endif

// Function prototypes
void SaveCanvases(TObjArray* canvases, const char* fileName);
void SaveCanvasesFromFile(const char* rootFile, 
			  const char* targetDir, 
			  const char* tag, 
			  const char* fileType);
TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir="");
int PrintPDFs(TObjArray* cList, TString dir="./", TString opt="");
int PrintPDF(TObjArray* cList, TString base, TString opt="pdf");
TCanvas* DrawObject(TObject* obj, 
		    TString drawopt="", 
		    TString title="",
		    TObjArray* cList = 0,
		    double xpx=700, double ypx=500);
void SetHistProps(TH1* h,
		  Int_t linecolor = kBlack,
		  Int_t fillcolor = kNone,
		  Int_t markercolor = kBlack,
		  Int_t markerstyle = kDot,
		  Double_t markersize = 1.0); 
void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t fillcolor,
		   Int_t markercolor,
		   Int_t markerstyle=kFullCircle,
		   Double_t markersize=1.0);
TGraphTime* Animation(TObjArray* moveObjs, 
		      TObjArray* statObjs,
		      TString opt="",
		      int sleeptime=50);
TGraphTime* Animation(TH2* h, 
		      TObjArray* statObjs=0,
		      TString opt="", 
		      int sleeptime=50,
		      int color=kBlack, 
		      int mkr=kFullCircle, 
		      double mkrsize=1.0,
		      double x1=0,double y1=0,double x2=0,double y2=0);
TGraphTime* Animation(TH3* h, 
		      TString propt="yx", // project z onto x-y plane 
		      int sleeptime=50,
		      TString opt="colz");

TMultiGraph* MultiGraph(TH2* h, TString opt="");

// Function definitions
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

int PrintPDFs(TObjArray* cList, TString dir, TString opt)
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

    if (gSystem->FindFile(dir.Data(), fileName.Data()))
      gROOT->Info("PrintPDFs()", "Overwriting %s", fileName.Data());
    
    if (!dir.EndsWith("/"))
      dir.Append("/");
    
    fileName.Prepend(dir);
    
    if (0)
      Info("UtilFns - PrintPDFs()", 
	   "dir = %s, fileName = %s", dir.Data(), fileName.Data());
    
    c->Print(fileName.Data());
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
      c = (TCanvas*)obj;
    else {
      gROOT->Warning("PlotUtils::print_pdf()", 
		     "list contains non-canvas object");
      continue;
    }
    if (!c) {
      gROOT->Error("PlotUtils::print_pdf()", "no canvas!");
      return -1;
    }

    // Slide numbering
    if (opt.Contains("number") || opt.Contains("#")) {
      c->cd();
      ltx.DrawLatex(0.95, 0.01, Form("%d", nPrinted+1));
    }
    
    // Multipage ps and pdf files
    if (nPrinted==0) {
      psOut = base + ext + "[";
      c->Print(psOut.Data()); // opens ps but doesn't print it
    }
    if (i < cList->GetEntries()) {
      psOut = base + ext;
      if(ext.Contains("pdf")) {
	TString pageName = Form("Title:%s", c->GetName());
	c->Print(psOut.Data(), pageName.Data());
      }
      else
	c->Print(psOut.Data());
    }
    if (i==cList->GetEntries()-1) {
      psOut = base + ext + "]"; // closes ps but doesn't print it
      if(ext.Contains("pdf")) {
	TString pageName = Form("Title:%s", c->GetName());
	c->Print(psOut.Data(), pageName.Data());
      }
      else
	c->Print(psOut.Data());
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

  if (!title.IsNull()) {
    c->SetName(title.Data());
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

void CopyProps(TObject* obj, TObjArray* arr)
{
  int lc=0, fc=0, mc=0, ms=0, mz=0;
  TH1    *h=0;
  TGraph *g=0;

  bool isTH1 = obj->InheritsFrom(TH1::Class());
  bool isTGr = obj->InheritsFrom(TGraph::Class());

  if (isTH1) {
    h  = dynamic_cast<TH1*>(obj);
    lc = h->GetLineColor();
    fc = h->GetFillColor();  
    mc = h->GetMarkerColor();
    ms = h->GetMarkerStyle();
    mz = h->GetMarkerSize();
  }
  else if (isTGr) {
    g  = dynamic_cast<TGraph*>(obj);
    lc = g->GetLineColor();
    fc = g->GetFillColor();  
    mc = g->GetMarkerColor();
    ms = g->GetMarkerStyle();
    mz = g->GetMarkerSize();
  }
  else {
    Warning("UtilFns - CopyProps()",
  	    "Class %s not recognized", obj->ClassName());
    return;
  }
  
  for (int i=0; i<arr->GetEntries(); i++) {
    
    if ((arr->At(i))->InheritsFrom(TH1::Class())) {
      h = (TH1*)arr->At(i);
      h->SetLineColor(   lc );
      h->SetFillColor(   fc );
      h->SetMarkerColor( mc );
      h->SetMarkerStyle( ms );
      h->SetMarkerSize(  mz );
    }
    else if ((arr->At(i))->InheritsFrom(TGraph::Class())) {
      g = (TGraph*)arr->At(i);
      g->SetLineColor(   lc );
      g->SetFillColor(   fc );
      g->SetMarkerColor( mc );
      g->SetMarkerStyle( ms );
      g->SetMarkerSize(  mz );
    }
  }

}

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t fillcolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize) 
{
  g->SetLineColor(linecolor);
  g->SetFillColor(markercolor);
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

    TObject* mov = moveObjs->At(n);
    
    // Add changing objects
    anim->Add(mov, n, opt);
    
    // Add stationary objects
    if (statObjs) {
      for (int i=0; i<statObjs->GetEntries(); i++) {
	TObject* sta = statObjs->At(i);
	anim->Add(sta, n, opt+"same");
      }
    }
    
  }

  anim->SetSleepTime(sleeptime); // ms (default = 0)

  // Rename auto-generated frame histo
  TH1* hf = (TH1*)gDirectory->Get("frame");
  hf->SetName(Form("frame%d", iAnim));
  
  return anim;
}

TGraphTime* Animation(TH2* h, 
		      TObjArray* statObjs,
		      TString opt, 
		      int sleeptime, 
		      int color, 
		      int mkr, 
		      double mkrsize,
		      double x1,double y1,double x2,double y2)
{
  static int id=0; id++;
  TObjArray* a = new TObjArray();

  for (int k=1; k<=h->GetNbinsY(); k++) {
    TH1D* hp = h->ProjectionX(Form("h%d_%d",id,k),k,k);
    
    if (x2>x1) hp->GetXaxis()->SetRangeUser(x1,x2);
    if (y2>y1) hp->GetYaxis()->SetRangeUser(y1,y2);

    //    hp->SetTitle(Form("%s, #lambda = %.2g",hp->GetTitle(), h->GetYaxis()->GetBinCenter(k)));
    //    hp->SetTitle(Form("%s, #lambda bin %d",hp->GetTitle(), k));

    SetHistProps(hp,color,0,color,mkr,mkrsize);
    hp->SetLineWidth(2);
    a->Add(hp);
  }
  return Animation(a, statObjs, opt, sleeptime);
}

// TODO: parse propt to determine automatically which proj to do according to Project3D options
TGraphTime* Animation(TH3* h,
		      TString propt, 
		      int sleeptime, 
		      TString opt)
{
  static int id=0; id++;
  TObjArray* a = new TObjArray();
  
  TAxis* axis = h->GetZaxis();
  
  for (int k=1; k<=axis->GetNbins(); k++) {
    axis->SetRange(k,k);
    TH2D* h2 = h->Project3D(Form("h2_%d_%d_%s",id,k,propt.Data()));
    h2->SetTitle(Form("k = %d",k));
    a->Add(h2);
  }
  return Animation(a, 0, opt, sleeptime);
}

TMultiGraph* MultiGraph(TH2* h, TString opt)
{
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  TMultiGraph* mg = new TMultiGraph();
  TGraph* g = 0;
 
  for (int i=0; i<ny; i++) {
    g = new TGraph(nx);
    g->SetLineWidth(2);
    g->SetLineColor(kBlue);
    g->SetTitle(Form("g%d",i));
    for (int j=0; j<nx; j++) {
      double x = h->GetXaxis()->GetBinCenter(j+1);
      double y = h->GetBinContent(j+1,i+1); 
      g->SetPoint(j,x,y);
    }
    mg->Add(g);
  }

  return mg;
}

void MakeBeamerSlides(TString dir, TString texFileName)
{
  // Create Beamer slides from a collection of images found in dir.

  TString extensions[3] = {"pdf", "png", "jpg"};
  TSystemDirectory tsd(dir.Data(), dir.Data());
  TList* files = tsd.GetListOfFiles();

  gSystem->RedirectOutput(texFileName, "w", 0);

  if (files) {
    TSystemFile* file = 0;
    TString fileName;
    TIter next(files);
    while((file = (TSystemFile*)next())) {
      fileName = file->GetName();
      bool isImage = false;
      for (int i=0; i<3; i++) {
	if (fileName.EndsWith(extensions[i].Data()))
	  isImage = true;
      }
      if (!file->IsDirectory() && isImage) {
	PrintSlide(fileName);
      }
    }
  }
  
  gSystem->RedirectOutput(0); // Back to stdout, stderr
  return;
}

void PrintSlide(TString fig)
{
  Printf("\\begin{frame}{title}{subtitle}");
  Printf("\\begin{center}");
  Printf("\\includegraphics[width=0.8\\textwidth]{%s}",fig.Data());
  Printf("\\end{center}");
  Printf("\\end{frame}\n\n");
}
