// MatrixUtils.C
// Matrix functions to augment ROOT's capabilities.
//
// - Matrix and vector <--> histogram conversions
// - Convenience methods for column/row manipulation
// - QR, CS, and Generalized SV decompositions
// - A few matrix-like operations for 2D histograms
// - A couple commonly-used matrix types
//
// To use in a ROOT script, either compiled or interpreted:
//   #include "MatrixUtils.C"
//   Then qualify names globally:
//     using MatrixUtils
//   Or not:
//     MatrixUtils::Hist2Vec(h)
//
// Please report issues to Andrew Adare andrewadare@gmail.com

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TDecompQRH.h>
#include <TDecompChol.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>

using std::cout;
using std::endl;

namespace MatrixUtils
{

struct QRDecompResult
{
  TMatrixD Q;
  TMatrixD R;
};

struct CSDecompResult
{
  TMatrixD C;
  TMatrixD S;
  TMatrixD U;
  TMatrixD V;
  TMatrixD Z;
  TVectorD alpha;        // Diag(C). Size l (CSD) = n (GSVD).
  TVectorD beta;         // Diag(S). Size l (CSD) = n (GSVD).
};

struct GSVDecompResult   // Output from GSVD().
{
  // Result of quotient aka generalized SVD of A,B.
  // A = U*C*XT, B = V*S*XT.
  TMatrixD C;
  TMatrixD S;
  TMatrixD U;
  TMatrixD V;
  TMatrixD XT;

  TVectorD alpha;        // Diag(C). Size n. Nonincreasing.
  TVectorD beta;         // Diag(S). Size n. Nondecreasing.
  TVectorD gamma;        // alpha/beta. Nonincreasing.
};

// Conversion methods
TVectorD Hist2Vec(const TH1 *h, TString opt=""); // use "unc" to get error
TMatrixD Hist2Matrix(const TH2 *h);
TVectorD Graph2Vec(const TGraph *g);
TH1D *Vec2Hist(const TVectorD &v, Double_t x1, Double_t x2,
               TString name, TString title="");
TH2D *Matrix2Hist(TMatrixD &A, TString hName,
                  double x1, double x2, double y1, double y2);
TH2D *Matrix2Hist(TMatrixD &A, TString hName,
                  double xbins[], double ybins[]);
TH2D *Matrix2Hist(TMatrixD &A, TMatrixD &errA, TString hName,
                  double xbins[], double ybins[]);
TH2D *BandedDiagonalMatrix(TH1 *hDpt,
                           const int nMeas,
                           const int nTrue,
                           double xm1=0,
                           double xm2=0,
                           double xt1=0,
                           double xt2=0);

// Matrix decompositions
QRDecompResult QRDecomp(TMatrixD &A);
QRDecompResult QLDecomp(TMatrixD &A);
CSDecompResult CSDecomp(TMatrixD &Q1, TMatrixD &Q2);
CSDecompResult CSDecompQ1Taller(TMatrixD &Q1, TMatrixD &Q2);
GSVDecompResult GSVD(TMatrixD &A, TMatrixD &B);

// Utility methods
void ReverseColumns(TMatrixD &A); // Reverse all columns in A
void ReverseColumns(TMatrixD &A, int col1, int col2); // col1, col2 included
void ReverseRows(TMatrixD &A);
void ReverseVector(TVectorD &v);
void SwapColumns(TMatrixD &A, int col1, int col2);
void SwapElements(TVectorD &v, int j1, int j2);
void NormalizeXSum(TH2 *hA, TH1 *hN=0); // Modifies hA in-place
void NormalizeRows(TMatrixD &A, TVectorD &normto); // Modifies A in-place
TH2 *TH2Product(TH2 *hA, TH2 *hB, TString name);
TH2D *TH2Sub(TH2 *h, int bx1, int bx2, int by1, int by2, TString name);

TMatrixD MoorePenroseInverse(TMatrixD &A, double tol = 1e-15); // Uses SVD
TMatrixD Null(TMatrixD &A); // Columns form a basis for the null space of A
int Rank(TMatrixD &A); // Uses SVD
TMatrixD Toeplitz(int m1, int n1, double col[], double row[]); // Pass in first col & row
TMatrixD LMatrix(const int n, const int kind, double eps = 1.e-5); // Matrix to define smoothing seminorm
TMatrixD DerivativeMatrix(int n, int d);
TVectorD Ones(int n); // n-vector of 1's
TVectorD ElemMult(const TVectorD &x, const TVectorD &y); // Element-wise vector multiplication
TVectorD ElemDiv(const TVectorD &x, const TVectorD &y, double div0val = 0.);   // Element-wise vector division
TMatrixD MultRowsByVector(const TMatrixD &M, const TVectorD &v);
TMatrixD DivRowsByVector(const TMatrixD &M, const TVectorD &v, bool makeZeroIfNaN=true); // Divide M rows by v elementwise: R(i,j) = M(i,j) / v(j).
TMatrixD DivColsByVector(const TMatrixD &M, const TVectorD &v, bool makeZeroIfNaN=true); // Divide M columns by v elementwise: R(i,j) = M(i,j) / v(i).
TMatrixD OuterProduct(TVectorD a, TVectorD b); // a*b'
void LTSolve(TVectorD &result, const TMatrixD &L, const TVectorD &y);
void LSolve(TVectorD &result, const TMatrixD &L, const TVectorD &y,
            const TMatrixD &W, const TMatrixD &T);

TMatrixD
Hist2Matrix(const TH2 *h)
{
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  TMatrixD m(nx, ny);
  for (Int_t j=0; j<ny; j++)
  {
    for (Int_t i=0; i<nx; i++)
    {
      m(i,j) = h->GetBinContent(i+1,j+1);
    }
  }
  return m;
}

TH2D *
Matrix2Hist(TMatrixD &A, TString hName,
            double x1, double x2, double y1, double y2)
{
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D *h = new TH2D(hName.Data(),hName.Data(),m,x1,x2,n,y1,y2);

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      h->SetBinContent(i+1, j+1, A(i,j));
    }
  }

  return h;
}

TH2D *
Matrix2Hist(TMatrixD &A, TString hName,
            double xbins[], double ybins[])
{
  // xbins and ybins better have size m+1 and n+1, respectively
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D *h = new TH2D(hName.Data(),hName.Data(),m,xbins,n,ybins);

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      h->SetBinContent(i+1, j+1, A(i,j));
    }
  }

  return h;
}

TH2D *
Matrix2Hist(TMatrixD &A, TMatrixD &errA, TString hName,
            double xbins[], double ybins[])
{
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D *h = new TH2D(hName.Data(),hName.Data(),m,xbins,n,ybins);

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      h->SetBinContent(i+1, j+1, A(i,j));
      h->SetBinError(i+1, j+1, errA(i,j));
    }
  }

  return h;
}


TH1D *
Vec2Hist(const TVectorD &v, Double_t x1, Double_t x2,
         TString name, TString title)
{
  int nb = v.GetNoElements();
  TH1D *h = new TH1D(name.Data(), title.Data(), nb, x1, x2);

  for (int i=0; i<nb; i++)
  {
    h->SetBinContent(i+1, v(i));
  }
  return h;
}

TVectorD
Hist2Vec(const TH1 *h, TString opt)
{
  // Returns TVectorD of the bin contents of the input histogram
  int nb = h->GetNbinsX();
  TVectorD v(nb);
  if (!h) return v;
  double val = 0.;
  for (Int_t i= 0; i<nb; i++)
  {
    if (opt.Contains("unc"))
      val = h->GetBinError(i+1);
    else
      val = h->GetBinContent(i+1);
    v(i) = val;
  }
  return v;
}

TVectorD
Graph2Vec(const TGraph *g)
{
  // Return a TVectorD from a TGraph (or inherited classes)
  int nb = g->GetN();
  TVectorD v(nb);
  if (!g) return v;

  for (int i=0; i<nb; i++)
  {
    v(i) = g->GetY()[i];
  }
  return v;
}

void
LSolve(TVectorD &result, const TMatrixD &L, const TVectorD &y,
       const TMatrixD &W, const TMatrixD &T)
{
  // Computes  x = L_p*y
  // where L_p is the A-weighted generalized inverse of L.
  int p = L.GetNrows();
  int n = L.GetNcols();
  int ly = y.GetNoElements();

  if (ly != p)
    gROOT->Warning("LSolve()",
                   "Input vector length = %d != %d", ly, p);
  if (p==n)
  {
    TMatrixD Linv(TMatrixD::kInverted, L);
    TVectorD x = Linv * y;
    result = x;
    return;
  }

  if (p > n)
  {
    TMatrixD Ltmp = L;
    TMatrixD Linv = MoorePenroseInverse(Ltmp);
    TVectorD x = Linv * y;
    result = x;
    return;
  }
  else
  {
    TMatrixD L11(L); L11.ResizeTo(p,p);
    TMatrixD T11(T); T11.ResizeTo(n-p,p);
    TMatrixD L11inv = MoorePenroseInverse(L11);
    TVectorD xhat = L11inv * y;
    TVectorD xhat0 = xhat;
    xhat0.ResizeTo(n); // Append n-p 0's
    TVectorD x = xhat0 - W*T11*xhat;
    result.ResizeTo(x.GetNoElements());
    result = x;
  }
  return;
}

void
LTSolve(TVectorD &result, const TMatrixD &L, const TVectorD &y)
{
  // Computes  x = (L_p)'*y
  // where L_p is the generalized inverse of L (aka L_A^dagger)

  int p = L.GetNrows();
  int n = L.GetNcols();
  int ly = y.GetNoElements();

  result.Zero();

  if (p == n)
  {
    TMatrixD LT(L);
    LT.Invert();
    LT.T();
    result = LT * y;
    return;
  }
  if (p > n)
  {
    TMatrixD LT(L); LT.T();
    result = MoorePenroseInverse(LT) * y;
    return;
  }
  else
  {

    TVectorD y1(y);
    TMatrixD L11(L);
    y1.ResizeTo(p);
    L11.ResizeTo(p,p);

    // Calculate result
    L11.Invert();
    L11.T();
    result = L11 * y1;

    if (0)
    {
      Printf("p, n, ly = %d, %d, %d", p, n, ly);
      Printf("y, y1, result: %d, %d, %d",
             y.GetNoElements(),
             y1.GetNoElements(),
             result.GetNoElements());
      Printf("L11 (%d x %d)",
             L11.GetNrows(), L11.GetNcols());
    }

  }
  return;
}

TMatrixD
MoorePenroseInverse(TMatrixD &A, double tol)
{
  // Compute the Moore-Penrose (pseudo)inverse of A as
  // V*diag(1/sigma_i)*U' (Numerical Recipes eq 2.6.6)
  // Inverse singular values < tol are zeroed.
  int m = A.GetNrows(), n = A.GetNcols();
  int max = m;
  if (n>max)
    max = n;
  // Get SVD components
  TDecompSVD decomp;

  if (m<n)
  {
    TMatrixD B(A);
    B.ResizeTo(n,n);
    decomp.SetMatrix(B);
  }
  else
  {
    decomp.SetMatrix(A);
  }

  TMatrixD UT = decomp.GetU(); UT.T();
  TVectorD sig = decomp.GetSig();
  TMatrixD V = decomp.GetV();
  TMatrixD Siginv(max,max);

  V.ResizeTo(max,max);

  for (int i=0; i<n; i++)
  {
    if (sig(i) > tol)
      Siginv(i,i) = 1./sig(i);
    else
      Siginv(i,i) = 0.;
  }

  TMatrixD Ainv = (V*Siginv)*UT;
  Ainv.ResizeTo(n,m);

  return Ainv;
}

TVectorD
Ones(int n)
{
  TVectorD ones(n);
  for (int i=0; i<n; i++)
    ones(i) = 1.0;
  return ones;
}

TVectorD
ElemDiv(const TVectorD &x, const TVectorD &y, double div0val)
{
  int nx = x.GetNoElements();
  int ny = y.GetNoElements();
  if (nx != ny)
  {
    gROOT->Error("ElemDiv()", "mismatch nx=%d, ny=%d", nx, ny);
    gSystem->Exit(-1);
  }
  TVectorD result(nx);
  for (int i=0; i<nx; i++)
  {
    result(i) = (y(i) > 1e-15) ? x(i) / y(i) : div0val;
  }
  return result;
}

TVectorD
ElemMult(const TVectorD &x, const TVectorD &y)
{
  int nx = x.GetNoElements();
  int ny = y.GetNoElements();
  if (nx != ny)
  {
    gROOT->Error("ElemMult()", "mismatch nx=%d, ny=%d", nx, ny);
    gSystem->Exit(-1);
  }
  TVectorD result(nx);
  for (int i=0; i<nx; i++)
  {
    result(i) = x(i)*y(i);
  }
  return result;
}

TMatrixD
DivColsByVector(const TMatrixD &M, const TVectorD &v,
                bool makeZeroIfNaN)
{
  // Divide M columns by v elementwise: R(i,j) = M(i,j) / v(i).
  // I.e. each row i is scaled by 1/v(i).
  TMatrixD R(M);
  int m = R.GetNrows(), n = R.GetNcols();

  if (v.GetNoElements() != m)
    Error("DivColsByVector()",
          "nrows %d != vector size %d", m, v.GetNoElements());

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      if (v(i) != 0)
        R(i,j) /= v(i);
      else if (makeZeroIfNaN)
        R(i,j) = 0;
    }
  }
  return R;
}

TMatrixD
DivRowsByVector(const TMatrixD &M, const TVectorD &v,
                bool makeZeroIfNaN)
{
  // Divide M rows by v elementwise: R(i,j) = M(i,j) / v(j).
  // I.e. each column j is scaled by 1/v(j).
  TMatrixD R(M);
  int m = R.GetNrows(), n = R.GetNcols();

  if (v.GetNoElements() != n)
    Error("DivRowsByVector()",
          "nrows %d != vector size %d", n, v.GetNoElements());

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      if (v(j) != 0)
        R(i,j) /= v(j);
      else if (makeZeroIfNaN)
        R(i,j) = 0;
    }
  }
  return R;
}

TMatrixD
MultRowsByVector(const TMatrixD &M, const TVectorD &v)
{
  // Multiply M rows by v elementwise: R(i,j) = M(i,j) * v(j).
  // I.e. each column j is scaled by v(j).
  TMatrixD R(M);
  int m = R.GetNrows(), n = R.GetNcols();

  if (v.GetNoElements() != n)
    Error("MultRowsByVector()",
          "ncols %d != vector size %d", n, v.GetNoElements());

  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      R(i,j) *= v(j);
    }
  }
  return R;
}
TMatrixD
Toeplitz(int m1, int n1, double col[], double row[])
{
  // Return a Toeplitz matrix T constructed from col[] (size m+1) and
  // row[] (size n+1). T has dimensions (m+1) x (n+1). T is assigned
  // as
  //
  //                   c[0] = r[0], i==j (diagonal)
  //    t[i][j] =      c[i-j],      i>j  (lower left)
  //                   r[j-i],      i<j  (upper right)
  //
  // To create arrays to pass in, follow this pattern:
  // c[m+1] = {t[0], t[-1], t[-2], ..., t[-m]} (i.e. t[-i] = c[i])
  // r[n+1] = {t[0], t[1],  t[2],  ..., t[n]}  (i.e. t[j]  = r[j])

  if (col[0] != row[0])
  {
    Warning("Toeplitz()",
            "col[0] (%f) != row[0] (%f). Using col[0].",
            col[0], row[0]);
  }

  TMatrixD T(m1, n1);
  T(0, 0) = col[0];
  for (int i=0; i<m1; i++)
  {
    for (int j=0; j<n1; j++)
    {
      if (i==j)  T(i, j) = col[0];
      else if (i>j)  T(i, j) = col[i-j];
      else if (i<j)  T(i, j) = row[j-i];
    }
  }
  return T;
}

int
Rank(TMatrixD &A)
{
  // Not very efficient, could replace with a rank-revealing QR
  // algorithm.
  TMatrixD N = Null(A);
  int nullity = N.GetNcols();
  int rank = TMath::Min(A.GetNrows(), A.GetNcols()) - nullity;
  return rank;
}

TMatrixD
Null(TMatrixD &A)
{
  // Return matrix whose columns form the orthonormal basis for the
  // nullspace of A. Perform SVD on A so that A = USV'. The columns of
  // V whose corresponding singular values are "zero" are the null basis
  // vectors.
  // TODO: pass in tolerance instead of hardcoding 1e-16
  int m = A.GetNrows(), n = A.GetNcols();

  // Result
  TMatrixD W(0, 0);

  // Get SVD components
  TDecompSVD decomp;

  if (m<n)
  {
    TMatrixD B = A;
    B.ResizeTo(n,n);
    decomp.SetMatrix(B);
  }
  else
  {
    decomp.SetMatrix(A);
  }

  TVectorD sig = decomp.GetSig();
  TMatrixD V = decomp.GetV();

  if (0)
  {
    sig.Print();
    V.Print();
  }

  int ndim = 0;
  for (int i=0; i<sig.GetNoElements(); i++)
  {
    if (sig(i) < 1e-16)
    {
      ndim++;
      W.ResizeTo(n, ndim);
      TMatrixDColumn(W,ndim-1) = TMatrixDColumn(V, i);
    }
  }

  return W;
}

TH2 *
TH2Product(TH2 *hA, TH2 *hB, TString name)
{
  // Return the matrix product hA*hB as a new TH2
  int nA = hA->GetNbinsY(), mB = hB->GetNbinsX();

  if (nA != mB)
  {
    Error("TH2Product()",
          "Number of A cols (%d) != B rows (%d)", nA, mB);
    return 0;
  }

  double xa1 = hA->GetXaxis()->GetXmin(), xa2 = hA->GetXaxis()->GetXmax();
  double yb1 = hB->GetYaxis()->GetXmin(), yb2 = hB->GetYaxis()->GetXmax();

  TMatrixD A = Hist2Matrix(hA);
  TMatrixD B = Hist2Matrix(hB);
  TMatrixD AB(A*B);
  return Matrix2Hist(AB, name, xa1, xa2, yb1, yb2);
}

TMatrixD
OuterProduct(TVectorD a, TVectorD b)
{
  // Return a*b'
  int na = a.GetNrows();
  int nb = b.GetNrows();
  TMatrixD M(na,nb);
  for (int i=0; i<na; i++)
  {
    for (int j=0; j<nb; j++)
    {
      M(i,j) = a(i)*b(j);
    }
  }
  return M;
}

TH2D *
TH2Sub(TH2 *h, int bx1, int bx2, int by1, int by2, TString name)
{
  // Return a sub-range of h (bx1..bx2) x (by1..by2), bx2,by2 included.
  double xw = h->GetXaxis()->GetBinWidth(1);
  double yw = h->GetYaxis()->GetBinWidth(1);
  double x1 = (bx1-1)*xw, x2 = bx2*xw;
  double y1 = (by1-1)*yw, y2 = by2*yw;
  int nx = bx2-bx1+1, ny = by2-by1+1;

  if (bx1<1 || bx2>h->GetNbinsX())
  {
    Error("TH2Sub()",
          "Requested x bins %d,%d out of range (1,%d)", bx1,bx2,h->GetNbinsX());
    return 0;
  }
  if (by1<1 || by2>h->GetNbinsY())
  {
    Error("TH2Sub()",
          "Requested y bins %d,%d out of range (1,%d)", by1,by2,h->GetNbinsY());
    return 0;
  }

  TH2D *hSub = new TH2D(name.Data(),name.Data(),nx,x1,x2,ny,y1,y2);

  int ipx=1, ipy=1;
  for (int ix=bx1; ix<=bx2; ix++)
  {
    for (int iy=by1; iy<=by2; iy++)
    {
      double val = h->GetBinContent(ix, iy);
      hSub->SetBinContent(ipx, ipy, val);
      ipy++;
    }
    ipy=1;
    ipx++;
  }

  return hSub;
}

TMatrixD
DerivativeMatrix(int n, int d)
{
  // Return a d-th derivative operator matrix.
  // n is the number of columns (as usual).
  // Default size is n-d x n.
  // Unit matrix is returned if d=0.

  int nd = n-d;
  TMatrixD L(nd, n);
  TVectorD c(d+1); // Use to populate L
  c.Zero();

  // Assign c
  if (d==0)
    c(0) = 1;
  else
  {
    c(0) = -1;
    c(1) = 1;
  }
  for (int i=2; i<=d; i++)
  {
    TVectorD a(d+2), b(c);
    a.SetSub(1, c);
    a.ResizeTo(d+1);
    b(d) = 0;
    c = a-b;
  }

  for (int i=0; i<nd; i++)
  {
    for (int j=0; j<d+1; j++)
    {
      L(i,j+i) = c(j);
    }
  }

  return L;
}

void
NormalizeXSum(TH2 *hA, TH1 *hN)
{
  // Normalize x-rows of hA so each row sums to 1.0 (default), or
  // optionally to the value of the jth bin of hN.

  int nx = hA->GetNbinsX();
  int ny = hA->GetNbinsY();

  if (hN)
    if (ny != hN->GetNbinsX())
      Error("NormalizeXSum()",
            "ny=%d != %d in hN", nx, hN->GetNbinsX());

  // xsum(j) contains sum of x cells in row j
  TVectorD xsum(ny);
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
    {
      xsum(j) += hA->GetBinContent(i+1,j+1);
    }
  }

  // Change bin contents of hA to normalized value, which is 1.0 if hN
  // is not passed in.
  for (int j=0; j<ny; j++)
  {
    double a = hN ? hN->GetBinContent(j+1) : 1.;
    double w = (xsum(j) != 0.) ? a/xsum(j) : a;
    for (int i=0; i<nx; i++)
    {
      double val = w*hA->GetBinContent(i+1,j+1);
      hA->SetBinContent(i+1,j+1, val);
    }
  }
}

void
NormalizeRows(TMatrixD &A, TVectorD &normto)
{
  // Scale elements of row i so that they sum to to normto(i)
  for (int i=0; i<A.GetNrows(); i++)
  {
    TVectorD row = TMatrixDRow(A,i);
    TMatrixDRow(A,i) *= normto(i)/row.Sum();
  }
}

QRDecompResult
QRDecomp(TMatrixD &A)
{
  // Compute QR decomposition of A using Householder transformations
  QRDecompResult qr;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD R(A);
  TMatrixD Qj(m,m);

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++)
  {
    TVectorD col = TMatrixDColumn(R,j);
    TVectorD x = col.GetSub(j,m-1);

    int sign = (col(j)<0.)? -1. : 1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(0) += alpha;

    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;

    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(j,j,H);

    // Update Q and R
    Q = Q*Qj;
    R = Qj*R;
  }

  bool requirePositivePivotEntries = true;
  if (requirePositivePivotEntries)
  {
    int r = TMath::Min(m,n);
    for (int i=0; i<r; i++)
    {
      if (R(i,i)<0)
      {
        TMatrixDRow(R,i) *= -1;
        TMatrixDColumn(Q,i) *= -1;
      }
    }
  }

  // Store results
  qr.Q.ResizeTo(Q);
  qr.R.ResizeTo(R);
  qr.Q = Q;
  qr.R = R;

  return qr;
}

QRDecompResult
QLDecomp(TMatrixD &A)
{
  // Compute QL decomposition of A using Householder transformations
  QRDecompResult ql;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD L(A);
  TMatrixD Qj(m,m);

  // For sign manipulation
  TMatrixD U1(m,m);
  U1.UnitMatrix();
  TMatrixD L1(n,n);
  L1.UnitMatrix();

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++)
  {
    TVectorD col = TMatrixDColumn(L,n-j-1);
    TVectorD x = col.GetSub(0,nIter-j);

    int sign = (x(nIter-j)<0.)? -1. : 1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(nIter-j) += alpha;

    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;

    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(0,0,H);

    // Update Q and L
    Q = Q*Qj;
    L = Qj*L;
  }

  bool requirePositivePivotEntries = true;
  if (requirePositivePivotEntries)
  {
    int r = TMath::Min(m,n);
    int d = m-n;
    for (int i=0; i<r; i++)
    {
      int row=i,col=i;
      if (m<n)
        col -= d;
      else if (m>n)
        row += 1;

      if (L(row,col)<0)
      {
        TMatrixDRow(L,row) *= -1;
        TMatrixDColumn(Q,row) *= -1;
      }
    }
  }

  // Store results
  ql.Q.ResizeTo(Q);
  ql.R.ResizeTo(L);
  ql.Q = Q;
  ql.R = L;

  return ql;
}


CSDecompResult
CSDecomp(TMatrixD &Q1, TMatrixD &Q2)
{
  // ***************************************
  // Q1 shorter or equal to Q2 (m <= p case)
  // ***************************************
  bool debug = false;
  int m,p,l,q1,q2,r,n;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  r  = 0;

  if (m > p)
  {
    Error("CSDecomp()",
          "Q1 rows (%d) <= Q2 (%d) rows required.\nExiting.",m,p);
    gSystem->Exit(-1);
  }

  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult csd;
  TVectorD alpha(l);
  TVectorD beta(l);

  // 1.
  q1 = TMath::Min(m,l); // 3
  q2 = TMath::Min(p,l); // 4

  // 2. SVD of Q2: Q2 = VSZ'
  TDecompSVD svdQ2(Q2);
  TMatrixD V    = svdQ2.GetU();    // p x p
  TMatrixD Z    = svdQ2.GetV();    // l x l
  beta = svdQ2.GetSig();           // l

  // 3-5. Re-order V, S, Z cols
  ReverseColumns(V,0,q2-1);
  ReverseColumns(Z);
  ReverseVector(beta);

  // 6.
  for (int i=0; i<l-q2; i++)
  {
    alpha(i) = 1.0;
    beta(i)  = 0.0;
  }

  if (debug)
  {
    cout << "V: ";  V.Print();
    //    cout << "S: ";  S.Print();
    cout << "Z: ";  Z.Print();
    cout << "beta:"; beta.Print(); // non-decreasing
  }

  // 7.
  // Find r where beta(r) <= 1/sqrt(2) < beta(r+1).  C++ problem:
  // index != dimension! Need two variables, r and rdim, to resolve
  // the ambiguity when r=0. rdim is the number of betas below 0.707,
  // r is the index.
  double thr = 1./TMath::Sqrt(2.);
  int rdim = 0;
  r = -1;
  for (int i=0; i<m-1; i++)
  {
    if (beta(i) <= thr && beta(i+1) > thr)
    {
      r = i;
      break;
    }
  }
  if (r == -1)
  {
    r = 0;
  }
  else rdim = r+1;

  if (debug)
  {
    Printf("r = %d, rdim = %d.",r,rdim);
  }

  // 8.
  TMatrixD T = Q1*Z; // (m x l)

  // 9.
  // QR decomp of T: T = UR
  QRDecompResult qrT = QRDecomp(T);
  TMatrixD U = qrT.Q;
  TMatrixD R = qrT.R;

  if (debug)
  {
    cout << "T: ";  T.Print();
    cout << "U and R: ";
    U.Print();
    R.Print();
  }

  // Get R2 and R3 from R
  TMatrixD R2 = R.GetSub(l-q2,r,l-q2,r);
  TMatrixD R3 = R.GetSub(rdim,q1-1,rdim,l-1);
  int r3r = R3.GetNrows();
  int r3c = R3.GetNcols();
  if (r3r < r3c)
    R3.ResizeTo(r3c, r3c);

  if (debug)
  {
    cout << "R2: ";  R2.Print();
    cout << "R3: ";  R3.Print();
  }

  // 10.
  // Compute SVD of R3: R3 = Ur*Cr*Zr'
  TDecompSVD svd2(R3);
  TMatrixD Ur = svd2.GetU();
  TMatrixD Zr = svd2.GetV();
  TVectorD a3 = svd2.GetSig();

  for (int i=0; i<a3.GetNrows(); i++)
    alpha(rdim+i) = a3(i);

  // 11.
  for (int i=q1; i<l; i++)
  {
    alpha(i) = 0.0;
    beta(i)  = 1.0;
  }

  // 12.
  for (int i=l-q2; i<rdim; i++)
  {
    alpha(i) = R2(i,i);
  }

  if (debug)
  {
    cout << "alpha: ";  alpha.Print();
  }

  // 13.
  // Form final U matrix
  // First resize U to undo TDecompSVD-required modification
  if (r3r < r3c)
  {
    Ur.ResizeTo(r3r,r3r);
  }

  TMatrixD Ur1(U);
  Ur1.UnitMatrix();
  Ur1.SetSub(rdim,rdim,Ur);

  if (debug)
  {
    cout << "Ur: ";  Ur.Print();
    cout << "Ur1: ";  Ur1.Print();
  }

  // 14.
  // Form final Z matrix
  TMatrixD Zrt(Zr);
  Zrt.ResizeTo(q2-rdim,q2-rdim);

  TMatrixD Zr1(Z);
  Zr1.UnitMatrix();
  Zr1.SetSub(rdim,rdim,Zr);

  if (debug)
  {
    cout << "Zr1: ";  Zr1.Print();
  }

  U = U*Ur1;
  Z = Z*Zr1;

  // 15.
  TMatrixD St(Zrt);
  St.Zero();
  for (int i=0; i<St.GetNcols(); i++)
    St(i,i) = beta(i+rdim);

  TMatrixD W = St*Zrt;

  // 16.
  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;
  n = TMath::Min(rdim,l-q2);

  // 17.
  TMatrixD Vpost(V);
  Vpost.UnitMatrix();
  Vpost.SetSub(rdim-n, rdim-n, Qw);
  V = V*Vpost;

  // Construct C from alpha
  for (int i=0; i<q1; i++)
    C(i,i) = alpha(i);

  // And S from beta
  for (int i=0; i<q2; i++)
    S(i,i) = beta(i);

  if (debug)
  {
    cout << "St: ";  St.Print();
    cout << "Zrt: ";  Zrt.Print();
    cout << "W: ";  W.Print();
    cout << "Qw: ";  Qw.Print();
    cout << "Vpost: ";  Vpost.Print();

    cout << "C: ";  C.Print();
    cout << "S: ";  S.Print();
    TMatrixD CC(C,TMatrixD::kTransposeMult,C);
    TMatrixD SS(S,TMatrixD::kTransposeMult,S);
    TMatrixD one = CC + SS;
    cout << "C'C + S'S: ";  one.Print();

    TMatrixD G(S,TMatrixD::kMultTranspose,Z);
    TMatrixD Ginv = MoorePenroseInverse(G);
    TMatrixD Vtest = Q2*Ginv;
    cout << "Vtest: ";  Vtest.Print();
    cout << "V: ";  V.Print();

    TMatrixD zz = Q2 - Vtest*G;
    cout << "zz: ";  zz.Print();
  }

  csd.C.ResizeTo(C);
  csd.S.ResizeTo(S);
  csd.U.ResizeTo(U);
  csd.V.ResizeTo(V);
  csd.Z.ResizeTo(Z);
  csd.alpha.ResizeTo(l);
  csd.beta.ResizeTo(l);

  csd.C = C;
  csd.S = S;
  csd.U = U;
  csd.V = V;
  csd.Z = Z;
  csd.alpha = alpha;
  csd.beta = beta;

  if (debug)
    Printf("m=%d, p=%d, l=%d, q1=%d, q2=%d, r=%d, n=%d",  m,p,l,q1,q2,r,n);
  return csd;
}

CSDecompResult
CSDecompQ1Taller(TMatrixD &Q1, TMatrixD &Q2)
{
  // ***************************************
  // m > p case
  // ***************************************
  bool debug = false;
  bool verbose = false;
  int m,p,l,q1,q2,r,rdim;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  rdim = 0;  // number of alpha values < 0.707
  r  = -1;   // index of 1st alpha < 0.707
  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult cs;
  TVectorD alpha(l);
  TVectorD beta(l);

  // 1.
  q1 = TMath::Min(m,l);
  q2 = TMath::Min(p,l);

  // 2.
  // SVD of Q1: UCZ'
  if (verbose)
    Info("CSDecompQ1Taller()",
         "Computing SVD of Q1 (%d x %d)",m,l);
  TDecompSVD svdQ1(Q1);
  TMatrixD U     = svdQ1.GetU();    // m x m
  TMatrixD Z     = svdQ1.GetV();    // l x l
  alpha          = svdQ1.GetSig();  // l

  // 3.
  for (int i = 0; i<l; i++)
  {
    beta(i)  = 1.;
    if (i>=q1) alpha(i) = 0.;
  }

  if (verbose)
    Info("CSDecompQ1Taller()",
         "m=%d, p=%d, l=%d, q1=%d, q2=%d",  m,p,l,q1,q2);
  if (debug)
  {
    Printf("\nSVD: Q1 = U*diag(alpha)*Z\'");
    cout << "U: ";  U.Print();
    cout << "alpha: ";  alpha.Print();
    cout << "Z: ";  Z.Print();
  }

  // 4.
  // Find r where alpha(r) >= 1/sqrt(2) > alpha(r+1)
  // r is the index, rdim is the # of values.
  double thr = 1./TMath::Sqrt(2.);
  for (int i=0; i<l-1; i++)
  {
    // if (alpha(i) >= thr && alpha(i+1) < thr) {
    //   r = i;
    if (alpha(i) >= thr)
      r = i;
    if (alpha(i+1) < thr)
      break;
  }
  if (r == -1)
  {
    r = 0;
  }
  else rdim = r+1;

  if (debug)
    Printf("r = %d, rdim = %d",r,rdim);

  // 5.
  TMatrixD T = Q2*Z;

  // 6.
  // QL decomp of T: T = VL
  QRDecompResult vl = QLDecomp(T);

  TMatrixD V = vl.Q;
  TMatrixD L = vl.R;

  if (debug)
  {
    cout << "T = Q2*Z = V*L: ";  T.Print();
    cout << "V: ";  V.Print();
    cout << "L (before permutation): ";  L.Print();
    TMatrixD Tcheck(T);
    Tcheck -= V*L;
    Printf("T - V*L sum: %g",Tcheck.Sum());
  }

  // Create permutation matrix Pi; L = Pi*L.
  TMatrixD Pi(p,p);
  TMatrixD Iq2(q2,q2); Iq2.UnitMatrix();
  Pi.SetSub(0,p-q2,Iq2);
  if (p>q2)
  {
    int d = p-q2;
    TMatrixD Id(d,d); Id.UnitMatrix();
    Pi.SetSub(p-1,0,Id);
  }
  L = Pi*L;

  // And its inverse (for later)
  TMatrixD PiInv(TMatrixD::kInverted, Pi);

  if (debug)
  {
    cout << "Pi: ";  Pi.Print();
    cout << "PiInv: ";  PiInv.Print();
    cout << "Pi*L: ";  L.Print();
  }

  // Get [ L_11 L_12 ], call it L1
  //  int lastrow = (r>0)? r-1 : 0;
  int l1r = (r>0)? p-q2+r-1 : p-q2;
  int l1c = (r>0)? l-q2+r-1 : l-q2;

  if (l1r > L.GetNrows()-1)
    l1r = L.GetNrows()-1;
  if (l1c > L.GetNcols()-1)
    l1c = L.GetNcols()-1;

  if (debug)
    Printf("Getting L1 (%d x %d) from L (%d x %d)",
           l1r,l1c,L.GetNrows(),L.GetNcols());

  TMatrixD L1 = L.GetSub(p-q2, l1r, p-q2, l1c);

  // Get L_23 := L2
  int row1 = p-q2+r;
  if (row1 > L.GetNrows()-1)
    row1 = L.GetNrows()-1;
  if (row1 < 0) row1 = 0;
  int col1 = l-q2+r;
  if (col1 > L.GetNcols()-1)
    col1 = L.GetNcols()-1;
  if (col1 < 0) col1 = 0;

  if (L1.GetNrows() < L1.GetNcols())   // So TDecompSVD works
  {
    if (debug)
      Printf("Resizing L1 (%d x %d) --> (%d x %d)",
             L1.GetNrows(), L1.GetNcols(), L1.GetNcols(), L1.GetNcols());
    L1.ResizeTo(L1.GetNcols(), L1.GetNcols());
  }

  if (debug)
  {
    cout << "L1 = [L_11 L_12] = Vl*Sl*Zl\': ";  L1.Print();
    Printf("Getting L2... row %d to %d, col %d to %d",
           row1, L.GetNrows()-1, col1, L.GetNcols()-1);
  }

  TMatrixD L2 = L.GetSub(row1, L.GetNrows()-1, col1, L.GetNcols()-1);

  if (debug)
  {
    cout << "L2: ";  L2.Print();
  }

  // 7.
  if (verbose)
    Info("CSDecompQ1Taller()",
         "Computing SVD of L1 (%d x %d)",L1.GetNrows(),L1.GetNcols());
  TDecompSVD svdl(L1);
  TMatrixD Vlbig = svdl.GetU();

  int rmax = TMath::Min(r-1, V.GetNrows()-1);

  if (debug)
    Printf("Getting Vl from Vlbig (%d x %d)... row 0 to %d, col 0 to %d",
           Vlbig.GetNrows(), Vlbig.GetNcols(), rmax,rmax);

  TMatrixD Vl = Vlbig.GetSub(0,rmax,0,rmax);
  TVectorD bl = svdl.GetSig();
  TMatrixD Zl = svdl.GetV();

  // 8-10.
  ReverseVector(bl);
  ReverseColumns(Vl);
  ReverseColumns(Zl);

  // 11.
  int imax = TMath::Min(l, bl.GetNrows());
  for (int i=0; i<imax; i++)
    beta(i) = bl(i);
  int nl2 = L2.GetNrows();
  for (int i=0; i<nl2; i++)
  {
    if (bl.GetNrows()+i < l)
      beta(bl.GetNrows()+i) = L2(i,i);
  }
  for (int i=0; i<l-q2; i++)
  {
    alpha(i) = 1.;
    beta(i)  = 0.;
  }

  // 12. Create V~ to post-multiply V
  TMatrixD Vpost(V);
  Vpost.UnitMatrix();
  Vpost.SetSub(p-q2,p-q2,Vl);

  // 13.
  if (debug)
    Printf("V = V (%d x %d) * Vpost( %d x %d) * PiInv(%d x %d)",
           V.GetNrows(), V.GetNcols(),
           Vpost.GetNrows(), Vpost.GetNcols(),
           PiInv.GetNrows(), PiInv.GetNcols());
  V = V*Vpost*PiInv;

  if (debug)
  {
    cout << "bl (after reversing):";  bl.Print();
    cout << "beta:";  beta.Print();
    cout << "Vl: ";  Vl.Print();
    cout << "Vpost: ";  Vpost.Print();
    cout << "V (= V*Vpost, final): ";  V.Print();
    cout << "Zl: ";  Zl.Print();
    cout << "Z: ";  Z.Print();
  }

  // 14.
  TMatrixD Zpost(Z);
  Zpost.UnitMatrix();
  Zpost.SetSub(0,0,Zl);
  Z = Z*Zpost;

  if (debug)
  {
    cout << "Z (final): "; Z.Print();
  }

  // 15. W = S~ * Zl
  TMatrixD W(Zl);
  W.Zero();
  for (int i=0; i<Zl.GetNcols(); i++)
    W(i,i) = alpha(i);

  if (debug)
  {
    cout << "S~: ";  W.Print();
  }

  W = W*Zl;
  if (debug)
  {
    cout << "W: ";  W.Print();
  }

  if (verbose)
    Info("CSDecompQ1Taller()",
         "Computing QRDecomp(W) (%d x %d)",W.GetNrows(),W.GetNcols());
  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;
  int ndiff =  U.GetNcols() - Qw.GetNrows();
  if (ndiff > 0)
  {
    int nr = Qw.GetNrows();
    Qw.ResizeTo(nr+ndiff, Qw.GetNcols()+ndiff);
    for (int j=0; j<ndiff; j++)
      Qw(nr+j, nr+j) = 1.;
  }
  if (debug)
  {
    cout << "Qw: ";  Qw.Print();
  }

  U = U*Qw;

  if (debug)
  {
    cout << "U: ";  U.Print();
  }

  // Finally, assign C and S matrices
  if (verbose)
    Info("CSDecompQ1Taller()","Formatting output...");
  for (int j=0; j<q1; j++)
    C(j,j) = alpha(j);
  for (int i=0; i<q2; i++)
    S(i,l-p+i) = beta(l-p+i);

  if (debug)
  {
    cout << "alpha: ";  alpha.Print();
    cout << "beta: ";   beta.Print();
    cout << "C: ";  C.Print();
    cout << "S: ";  S.Print();
    TMatrixD upper(C, TMatrixD::kMultTranspose, Z);
    TMatrixD lower(S, TMatrixD::kMultTranspose, Z);
    upper = U*upper;
    upper = Q1-upper;
    lower = V*lower;
    lower = Q2-lower;
    cout << "Q1 - U*C*Z\': ";  upper.Print();
    cout << "Q2 - V*S*Z\': ";  lower.Print();
  }

  cs.C.ResizeTo(C);
  cs.S.ResizeTo(S);
  cs.U.ResizeTo(U);
  cs.V.ResizeTo(V);
  cs.Z.ResizeTo(Z);
  cs.alpha.ResizeTo(l);
  cs.beta.ResizeTo(l);

  cs.C = C;
  cs.S = S;
  cs.U = U;
  cs.V = V;
  cs.Z = Z;
  cs.alpha = alpha;
  cs.beta = beta;

  return cs;
}


GSVDecompResult
GSVD(TMatrixD &A, TMatrixD &B)
{
  bool debug = false;
  GSVDecompResult g;
  int m,p,n,r;
  m = A.GetNrows();
  n = A.GetNcols();
  p = B.GetNrows();

  // TODO: check # A cols = B cols, and  m >= n >= p

  TMatrixD M(m+p,n);

  // 1. SVD of M = [A;B]: M = Q*S*Z'
  // TMatrixDSub(M, 0, m-1,   0, n-1) += A;
  // TMatrixDSub(M, m, m+p-1, 0, n-1) += B;
  M.SetSub(0, 0, A);
  M.SetSub(m, 0, B);

  if (M.GetNoElements() > 100000)
    Printf("GSVD(): "
           "Computing initial SVD on M (%d x %d)...", m+p, n);
  TDecompSVD svdM(M);
  TMatrixD Q  = svdM.GetU(); // m+p x m+p
  TMatrixD Z  = svdM.GetV(); //   n x n
  TVectorD sv = svdM.GetSig();

  r = sv.GetNrows(); // Assume M has full rank to save time
  //r = Rank(M);

  if (debug)
    Printf("Rank(M): %d. M = QSZ\', Q (%d x %d) Z (%d x %d)",
           r, Q.GetNrows(), Q.GetNcols(), Z.GetNrows(), Z.GetNcols());

  TMatrixD Sr(r,r); // Sigma_r submatrix from Sigma
  for (int i=0; i<r; i++)
    Sr(i,i) = sv(i);

  // 2. Partition Q to match dimensions of A and B.
  TMatrixD Q1 = Q.GetSub(0,m-1,0,n-1);
  TMatrixD Q2 = Q.GetSub(m, m+p-1, 0, n-1);
  bool q1Taller = Q1.GetNrows() > Q2.GetNrows();

  // 3. Do CS decomposition
  CSDecompResult csd =
    (q1Taller)? CSDecompQ1Taller(Q1, Q2) : CSDecomp(Q1, Q2);

  // 4. Assign output struct members
  g.U.ResizeTo(m,n); //csd.U);
  g.V.ResizeTo(csd.V);
  g.C.ResizeTo(n,n); //csd.C);
  g.S.ResizeTo(csd.S);
  g.alpha.ResizeTo(csd.alpha);
  g.beta.ResizeTo(csd.beta);
  g.gamma.ResizeTo(n);

  g.U = csd.U.GetSub(0,m-1,0,n-1);
  g.V = csd.V;
  g.C = csd.C.GetSub(0,n-1,0,n-1);
  g.S = csd.S;
  g.alpha = csd.alpha;
  g.beta = csd.beta;
  g.gamma = ElemDiv(g.alpha, g.beta, 9999e12);

  // Create X' from V'Sigma (upper left) and W (lower right)
  TMatrixD XT(n,n);
  XT.SetSub(0,0,TMatrixD(csd.Z,TMatrixD::kTransposeMult,Sr));
  for (int i=r; i<n; i++)
  {
    XT(i,i) = Sr(r-1,r-1);
  }
  XT = TMatrixD(XT, TMatrixD::kMultTranspose, Z);

  g.XT.ResizeTo(XT);
  g.XT = XT;

  return g;
}


void
ReverseColumns(TMatrixD &A, int col1, int col2)
{
  // Reverse the column sequence col1...col2, col1 and col2 included.
  TMatrixD B(A);
  int n = B.GetNcols();
  int ncols = col2-col1+1;
  if (ncols > n)
  {
    Error("ReverseColumns()",
          "Requested column range (%d-%d) > out of bounds (0-%d).",
          col1,col2,A.GetNcols());
    return;
  }

  for (int j=col1; j<=col2; j++)   // j is col index of B
  {
    TMatrixDColumn(B, j) = TMatrixDColumn(A, ncols-j-1);
  }
  A = B;
}

void
ReverseColumns(TMatrixD &A)
{
  TMatrixD B(A);
  int n = B.GetNcols();

  for (int j=0; j<n; j++)
  {
    TMatrixDColumn(B, j) = TMatrixDColumn(A, n-j-1);
  }
  A = B;
}

void
ReverseRows(TMatrixD &A)
{
  TMatrixD B(A);
  int m = B.GetNrows();

  for (int i=0; i<m; i++)
  {
    TMatrixDColumn(B, i) = TMatrixDColumn(A, m-i-1);
  }
  A = B;
}

void
ReverseVector(TVectorD &v)
{
  TVectorD vtmp(v);
  int m = v.GetNrows();
  for (int i=0; i<m; i++)
  {
    vtmp(i) = v(m-i-1);
  }
  v = vtmp;
}

void
SwapColumns(TMatrixD &A, int col1, int col2)
{
  int nc = A.GetNcols();
  if (col1 >= nc || col2 >= nc)
    Error("SwapColumns", "col1 or col2 index out of bounds");

  TMatrixD B(A);

  TMatrixDColumn(B, col1) = TMatrixDColumn(A, col2);
  TMatrixDColumn(B, col2) = TMatrixDColumn(A, col1);

  A = B;
}

void
SwapElements(TVectorD &v, int j1, int j2)
{
  int nr = v.GetNrows();
  if (j1 >= nr || j2 >= nr)
    Error("SwapElements", "an index is out of bounds");

  TVectorD v2(v);

  v2(j1) = v(j2);
  v2(j2) = v(j1);

  v = v2;
}

} // End namespace MatrixUtils