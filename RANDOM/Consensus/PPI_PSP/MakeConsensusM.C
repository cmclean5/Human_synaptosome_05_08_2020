//---
// 1)  module load root
// 2)  > root -l
// 3)  root> .L MakeConsensusM.C
// 4)  root> MakeConsensusM("full/path/to/files/",1,100)
//---

//---ROOT LIBRARIES
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrix.h"
#include "TMatrixD.h"

//---STD C++ LIBRARIES
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

inline void PlotStyle()
  {
    gROOT->SetStyle("Plain");
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetTitleColor(10);
    gStyle->SetOptStat("");
    gStyle->SetOptFit();
    gStyle->SetStatW(0.20);
    gStyle->SetStatH(0.13);
    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.46);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetStripDecimals(false);
    gStyle->SetHistTopMargin(0.2); 
    gStyle->SetPalette(1);

    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
  }

inline void print_message(const char *message){

  cout << TString(message) << endl;

}

// Returns the probability of x, given the Gaussian distribution described by mu and sigma.
double pdf(double x, double mu, double sigma)
{

  double pi = 3.14159265; 
  return exp( -1.0 * (x - mu) * (x - mu) / (2.0 * sigma * sigma)) / (sigma * sqrt(2.0 * pi));

}


// Returns the integral from -inf to x of any function that accepts x and 2 other parameters
void cdf(double x, double arg1, double arg2, double &res){
  
  double _sum  = 0.0;
  double ninf = -0.001;
  int    n    = 50;
  
  res         = 0.0;
  double res2 = 0.0;
  
  for (int k = 1; k < n-1; k++)
    _sum += pdf( x + (double)k*(x-ninf)/n , arg1, arg2);
  
  res = (double)(x - ninf/n) * ( pdf(x,arg1,arg2) + pdf(ninf,arg1,arg2)/2.0 + _sum );
  
}

void AreaUnderCurve(TMatrixD Ch, int K, TVectorD &auch, TVectorD &auc_tallyh){

  //--- The consensus Matrix
  int N    = Ch.GetNrows();
  int ndiv = 10;

  double _auc  = 0.0;
  double sum   = 0.0;
  double step  = 1.0/ndiv;

  double Xo    = 0.0;
  double Xi    = step;

  double mu    = 0.0;
  double sigma = 1.0;

  double norm  = N * ( N - 1.0 )/2.0;

  double _res  = 0.0;

  _auc         = 0.0;

  for(int m=2; m<= ndiv; m++){

    Xo  = Xi;
    Xi += step; 

    sum = 0.0;

    for(int i=0; i<N; i++){
      for(int j=i+1; j<N; j++){
	if( Ch(i,j) < Xi ){
	  _res = 0.0;
	  cdf( Xi, mu, sigma, _res );
	  sum += _res;
	}
      }
    }

    sum   = sum/norm;    
    _auc += fabs( Xi - Xo ) * sum;

  }

  auch[K]       += _auc;
  auc_tallyh[K] += 1.0;

  cout << auch[K] << " -- " << auc_tallyh[K] << " -- " << K << endl;

}

void calculateConsensusM(vector<int> _ki, vector<int> _kj, vector<int> _kk, TMatrixD &Ih, TMatrixD &Mh, TMatrixD &Ch){

  //Calculate the Concensus matrix
  int NG = (int)_ki.size();
  TMatrixD Ii(NG,NG);
  TMatrixD Mi(NG,NG);
  for(int i=0; i<_ki.size(); i++){
	    
    int comi = _kk[i];
    int keyi = _kj[i];
    
    for(int j=i; j<_kj.size(); j++){
      
      int comj = _kk[j];
      int keyj = _kj[j];

      Ii(i,j) = 0;
      Mi(i,j) = 0;
      if( keyi != -1 && keyj != -1 ){
	if(keyi == keyj)
	  Ii(i,j) = 0.5;
	else 
	  Ii(i,j) = 1;
		
	if(comi == comj && comi != -1 && comj != -1){
	  if( keyi == keyj )
	    Mi(i,j) = 0.5;
	  else
	    Mi(i,j) = 1;
	}
      }
	    
      
    }
  }

  Ih  = Ii;
  Ih += Ii.T();
	  
  Mh  = Mi;
  Mh += Mi.T();

  for(int i=0; i<NG; i++){
    for(int j=0; j<NG; j++){
      if(Ih(i,j) == 0){
	Ch(i,j) = 0;
      }else {
	Ch(i,j) = Mh(i,j)/Ih(i,j);
      }
    }
  }

}


int indicatorFunc(int ele_i, TVectorD Ik){

  for(int i=0; i<Ik.GetNoElements(); i++){
    if(Ik(i) == ele_i)
      return 1;
  }

  return 0;

}

bool storeElement( TMatrixD &mat_c, TMatrixD &mat_ids, int j, int ids, int c){

  if( c == -1){
    //for(int i=0; i<mat_c.GetNcols(); i++){
    //if(mat_c(0,i) == j)
    //return true;
      
    //if( mat_c(0,i) == 0){
    //mat_c(0,i) = j;
    //return true;
    //}
    //}
    return true;
  }

  if( c < mat_c.GetNrows() && c != -1 ){
    for(int i=0; i<mat_c.GetNcols(); i++){
      if(mat_c(c,i) == j)
	return true;

      if(mat_c(c,i) == 0){
	mat_c(c,i)   = j;
	mat_ids(c,i) = ids;
	return true;
      }
    }
  }

  if( c > mat_c.GetNrows() ){
    //cout << "c" << c << endl;
    mat_c.ResizeTo((c+1),500);
    mat_ids.ResizeTo((c+1),500);
    mat_c(c,0)   = j;
    mat_ids(c,0) = ids; 
    return true;
  }

  return true;

}

void MakeConsensusM(const char *path, int startJobID, int endJobID, const char *dataFile) {

  TMatrixD clusters(30,500);
  TMatrixD clusterIDs(30,500);
  clusters.Zero();
  clusterIDs.Zero();

  vector<int> key_listi;
  vector<int> key_listj;
  vector<int> key_listk;
  
  bool initM = true;
  TMatrixD I, M, C;

  int NJobs = 0;  
  int tally = startJobID;
  int stop  = 1000;

  int max_com       = 0;
  int min_com       = 100;
  int max_com_i     = 0;
  int max_com_tally = 0;

  TVectorD auc(500);
  TVectorD auc_tally(500);
  auc.Zero();
  auc_tally.Zero();
  
  char FileName[200];
  //sprintf(FileName,"consensusmatrix_%s.txt",dataFile);
  //fstream *fileout  = new fstream(FileName,ios_base::out);
  
  fstream *fileout  = new fstream("consensusmatrix.txt",ios_base::out);
  fstream *fileout2 = new fstream("auc.txt",ios_base::out);
  
  while( tally <= endJobID && NJobs < stop ){

    //char FileName[200];
    sprintf(FileName,"%s%d/consensusout.txt",path,tally);

    ifstream filein;
    filein.open(FileName);
    if( filein.is_open() ){

      char comments[256];
      filein.getline(comments, 256);
      //print_message(comments);

      int i,j,k;
      key_listi.clear();
      key_listj.clear();
      key_listk.clear();

      int ind   = 0;
      max_com_i = 0;
      //read data
      while (filein >> i >> j >> k){
	
	//--test 
	if( k == -1 ){
	  j = -1;
	} else {
	  //j = i;
	  i = j;
	}

	key_listi.push_back(i);
	key_listj.push_back(j);
	key_listk.push_back(k);

	if(k > max_com_i)           max_com_i = k; 
	if(k > max_com)             max_com   = k;
	if(k < min_com  && k != -1) min_com   = k;
	storeElement(clusters, clusterIDs, ind++, j, k);
      }
      
      max_com_tally += max_com_i;

      int N = key_listi.size();
      if(initM){
	cout << "N " << N << endl;
	I.ResizeTo(N,N);
	M.ResizeTo(N,N);
	initM = false;
      }

      //if( auc.GetNrows() == N ){
      //auc.ResizeTo(N);
      //auc_tally.ResizeTo(N);
      //}

      if( N == I.GetNrows() ){
	TMatrixD tempI(N,N);
	TMatrixD tempM(N,N);
	TMatrixD tempC(N,N);
	calculateConsensusM(key_listi, key_listj, key_listk, tempI, tempM, tempC);
	//AreaUnderCurve(tempC, max_com_i, auc, auc_tally);//-- comment out if searching for optimum K

	I += tempI;
	M += tempM;      

	NJobs++;
      }
      //cout << "Job completed" << endl;
      //cout << FileName << endl;

    } else {
      ;
      //cout << "Error opening file: " << endl;
      //cout << FileName << endl;

    }

    tally++;
    
  }

  //--- The consensus Matrix
  int N = I.GetNrows();
  C.ResizeTo(N,N);

  cout << "Matrix C" << endl;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(I(i,j) == 0){
	C(i,j) = 0;
      }else {
	C(i,j) = M(i,j)/I(i,j);
      }
    }
  }

  
  

  cout << "Number of good Jobs: " << NJobs << endl;
  cout << "max_com: " << max_com << " min_com: " << min_com << endl; 
  cout << "Avg. com: " << max_com_tally/NJobs << endl;

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(j+1 == N)
	(*fileout) << C(i,j) << endl;
      else
	(*fileout) << C(i,j) << ",";
    }
  }


  //cout << "printing auc V. community values " << endl;
  //for(int i=0; i<N; i++){
  //if( auc_tally[i] != 0 )
  //(*fileout2) << auc[i]/auc_tally[i] << "\t" << i+1 << endl;
  //else 
  //(*fileout2) << 0 << "\t" << i+1 << endl;
  //}

}
