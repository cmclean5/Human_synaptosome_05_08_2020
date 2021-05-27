#include "NetworkEnrichment.h"


NetworkEnrichment::NetworkEnrichment() : buildSets() {

  int i;

  this->freedMemory.push_back(false);

  this->baseNAME.push_back("permute_p_values");
  this->baseNAME.push_back("overlapCOMS");
  this->baseNAME.push_back("overlapNTWK");
  this->baseNAME.push_back("overlap");
  
  this->comSIZE   = 0;
  this->geneCOM   = 0;
  this->annoSIZE  = 0;
  
  this->studies   = 0;
  this->overlap   = 0;
  this->muab      = 0;
  this->muCab     = 0;
  this->nab       = 0;

  this->p_values  = 0;
  this->padjusted = 0;
  this->permute   = 0;  

  this->p_valuesT = 0;
  this->padjustedT= 0;
  this->permuteT  = 0;  
  
  this->p_valuesD = 0;
  this->padjustedD= 0;
  this->permuteD  = 0;

  this->p_valuesDT = 0;
  this->padjustedDT= 0;
  this->permuteDT  = 0;  

  this->p_dist     = 0;
  this->padjustedRD= 0;
  this->reldist    = 0;
  
  this->dels[0] = '\t';
  
  this->ANNOindex = 0;
  this->NoP       = 1;
  this->FDRtest   = 1;
  this->KOFFSET   = 0;
  this->isOFFSET  = false;
  this->printID   = false;
  this->printAn   = false;

  this->printMeanMu   = true;
  this->printFC       = false;
  this->printOneSided = false;
  this->printTwoSided = true;
  this->printALT      = true;
  this->calRelDist    = false;
  
  this->sigma[0] = 0.05;
  this->sigma[1] = 0.01;
  this->sigma[2] = 0.001;

  this->FDR       = 0.0;
  this->PV        = 0.0;
  this->SIGMA     = sigma[0];
  this->LEVEL     = 0;
  this->TESTS     = 0;

  this->pesudocount = 1.0; //to avoid p-values of zero, in the permutation test
  
  this->FDRmethods.push_back("BY");
  this->FDRmethods.push_back("BH");
  this->FDRmethods.push_back("BL");
  
  for( i=0; i<SIGMASIZE; i++ ){
    this->bonferroni[i]=0;
  }

  for( i=0; i<OVERLAPSIZE; i++ ){
    this->MINOVERLAP[i]=0;
  }

  //---set gsl random number generator, using taus2 generator
  g = gsl_rng_alloc(gsl_rng_taus2);
  seed = 0;
  seedOffset = 0;
  
}

NetworkEnrichment::NetworkEnrichment( vector<string> INfiles ) : buildSets() { 

  int i;
  int infilesSIZE=INfiles.size();
  
  const char* files[infilesSIZE];
  char        buffer[BUFFERSIZE];
  string      tmpdir;

  this->freedMemory.push_back(false);

  this->baseNAME.push_back("permute_p_values");
  this->baseNAME.push_back("overlapCOMS");
  this->baseNAME.push_back("overlapNTWK");
  this->baseNAME.push_back("overlap");
  
  this->comSIZE  = 0;
  this->geneCOM  = 0;
  this->annoSIZE = 0;
  
  this->studies   = 0;
  this->overlap   = 0;
  this->muab      = 0;
  this->muCab     = 0;
  this->nab       = 0;

  this->p_values  = 0;
  this->padjusted = 0;
  this->permute   = 0;  

  this->p_valuesT = 0;
  this->padjustedT= 0;
  this->permuteT  = 0;  
  
  this->p_valuesD = 0;
  this->padjustedD= 0;
  this->permuteD  = 0;

  this->p_valuesDT = 0;
  this->padjustedDT= 0;
  this->permuteDT  = 0; 

  this->p_dist     = 0;
  this->padjustedRD= 0;
  this->reldist    = 0;
  
  this->dels[0] = '\t';
  
  this->ANNOindex = 0;
  this->NoP       = 1000;
  this->FDRtest   = 1;
  this->KOFFSET   = 0;
  this->isOFFSET  = false;
  this->printID   = false;
  this->printAn   = false;

  this->printMeanMu   = true;
  this->printFC       = false;
  this->printOneSided = false;
  this->printTwoSided = true;
  this->printALT      = true;
  this->calRelDist    = false;
  
  this->sigma[0] = 0.05;
  this->sigma[1] = 0.01;
  this->sigma[2] = 0.001;

  this->FDR       = 0.0;
  this->PV        = 0.0;
  this->SIGMA     = sigma[0];
  this->LEVEL     = 0;
  this->TESTS     = 0;

  this->pesudocount = 1.0; //to avoid p-values of zero, in the permutation test
  
  this->FDRmethods.push_back("BY");
  this->FDRmethods.push_back("BH");
  this->FDRmethods.push_back("BL");
  
  for( i=0; i<SIGMASIZE; i++ )  { this->bonferroni[i]=0; }

  for( i=0; i<OVERLAPSIZE; i++ ){ this->MINOVERLAP[i]=0; }
   
  for(i=0; i<infilesSIZE; i++)  { files[i] = INfiles[i].c_str(); }
  
  vector<string>().swap(INfiles);//free space of vector

  //build our annotation sets
  addSets(files,infilesSIZE);  

  //printList( Clist, Clines);
    
  //Set default values for... but will reset/check this in each public calculation function
  //N -> Number of genes, network size
  //M -> Number of communities
  //F -> Number of annotation types
  N = Clines;
  M = COMS.size();
  F = Fsize[ANNOindex];
  
  //---set gsl random number generator, using taus2 generator
  g = gsl_rng_alloc(gsl_rng_taus2);
  seed = 0;
  seedOffset = 0;
  
}

NetworkEnrichment::~NetworkEnrichment(){

  freeMemory();
  
  gsl_rng_free(g);

}

//--------------------------------------------------------
//  BASE UTILITY FUNCTIONS
//--------------------------------------------------------
void NetworkEnrichment::freeMemory(){

  if( freedMemory[0] == false ){
  
  if(studies!=0) { free(studies);  }
  if(overlap!=0) { free(overlap);  }
  if(muab!=0)    { free(muab);     }
  if(muCab!=0)   { free(muCab);    }
  if(nab!=0)     { free(nab);      }

  if(p_values!=0)  { free(p_values);  }
  if(padjusted!=0) { free(padjusted); }
  if(permute!=0)   { free(permute);   }

  if(p_valuesT!=0) { free(p_valuesT); }
  if(padjustedT!=0){ free(padjustedT);}
  if(permuteT!=0)  { free(permuteT);  }
  
  if(p_valuesD!=0) { free(p_valuesD); }
  if(padjustedD!=0){ free(padjustedD);}
  if(permuteD!=0)  { free(permuteD);  }

  if(p_valuesDT!=0) { free(p_valuesDT); }
  if(padjustedDT!=0){ free(padjustedDT);}
  if(permuteDT!=0)  { free(permuteDT);  }

  if(p_dist!=0)      { free(p_dist); }
  if(padjustedRD!=0) { free(padjustedRD); }
  if(reldist!=0)     { free(reldist); }
    
  
  if( comSIZE!=0 ) { free(comSIZE);  }
  if( geneCOM!=0 ) { free(geneCOM);  }
  if( annoSIZE!=0 ){ free(annoSIZE); }

  freedMemory[0] = true;
  
  }

}

int NetworkEnrichment::getKOffset (){ return KOFFSET; }

void NetworkEnrichment::setKOffset ( int offset ){

  if( Mmin == 1 ){ return; }

  //We'll only use if Mmin is '0',
  //and keep offset set to '1'. 
  if( Mmin == 0 ){
  
    if( offset < 0 ){ offset = 0; }
    
    KOFFSET  = offset;
    isOFFSET = true;

    addKOffset( KOFFSET );

    M = COMS.size();

  }
  
    
}

void NetworkEnrichment::setNoP ( int perms ){

  if( (perms > 0) && (perms <= this->NoP) ){
    this->NoP = perms;
  }
  
}


void NetworkEnrichment::setPesudoCount ( double offset ){

  if( (offset >= 0) && (offset <= this->NoP) ){
    this->pesudocount = offset;
  }
  
}

void NetworkEnrichment::setExpectedOverlap ( bool setMeanMu ){

  this->printMeanMu = setMeanMu;
  
}

void NetworkEnrichment::setFoldChange ( bool setFC ){

  this->printFC = setFC;
  
}

void NetworkEnrichment::setRelDist ( bool setRD ){

  this->calRelDist = setRD;
  
}

void NetworkEnrichment::oneSided (){

  this->printTwoSided = false;
  this->printOneSided = true;
  
}

void NetworkEnrichment::twoSided (){

  this->printTwoSided = true;
  this->printOneSided = false;
  
}

void NetworkEnrichment::setALT ( bool setAlt ){

  this->printALT = setAlt;
  
}


void NetworkEnrichment::setPrintID ( bool setPrintId ){

  this->printID = setPrintId;

  //cout << "> printID = " << this->printID << endl;
  
}

void NetworkEnrichment::setPrintAn ( bool setPrintAn ){

  this->printAn = setPrintAn;

  //cout << "> printID = " << this->printID << endl;
  
}


double NetworkEnrichment::getMINOVERLAP ( int INDEX ){

  if( INDEX >= 0 && INDEX < OVERLAPSIZE){
    return MINOVERLAP[INDEX];
    }
  
  return -1;
  
}


void NetworkEnrichment::setMINOVERLAP ( int INDEX, int VALUE ){

  if( INDEX >= 0 && INDEX < OVERLAPSIZE){
    if( VALUE >= 0 ){
      MINOVERLAP[INDEX] = VALUE;
    }
  }
  
}

void NetworkEnrichment::setANNOindex ( int INDEX ){

  if( INDEX >= 0 && INDEX < Alist.size() ){
    ANNOindex = INDEX;
    F         = Fsize[ANNOindex];
    cout << "ANNOindex: " << ANNOindex << endl;//" => " << files[(ANNOindex+1)] << endl;
  }
  
}

void NetworkEnrichment::seedOffSet (bool useSeed, int newSeed ){

  if( useSeed==true ){
    cout << "" << endl;
    cout << "> newSeed: " << newSeed << endl;
    seedOffset = (unsigned long int) newSeed;
  }


}
  

//---Initialize random seed:
void NetworkEnrichment::setSeed (bool useSeed, int newSeed ){

  unsigned long int x,y,z,c,t;

  c = (unsigned long int) 6543217;

  x = (unsigned long int) time(NULL);

  seedOffset != 0 ? y = (unsigned long int) seedOffset : y = (unsigned long int) 987654321;

  z = (unsigned long int) getpid();

  cout << "" << endl;
  cout << "----" << endl;
  cout << "setting random number seed:" << endl;
  cout << "> x = " << x << endl;
  cout << "> y = " << y << endl;
  cout << "> z = " << z << endl;
  cout << "----" << endl;

  //--- JKISS RGN
  x  = 314527869 * x + 1234567;
  y ^= y << 5; y ^= y >> 7; y ^= y << 22;
  t  = 4294584393 * z + c; c = t >> 32; z = t; 

  seed = (unsigned long int) (x + y + z);
  
  cout << "> setSeed: " << seed << endl;
  cout << "----" << endl;

  gsl_rng_set(g , seed);
  
}

void  NetworkEnrichment::setFDRmethod ( const char* Method ){

  if( strcmp( Method, FDRmethods[0].c_str() ) == 0 ){  FDRtest = 1; cout << "set BY FDR Method" << endl; }

  if( strcmp( Method, FDRmethods[1].c_str() ) == 0 ){  FDRtest = 2; cout << "set BH FDR Method" << endl; }

  if( strcmp( Method, FDRmethods[2].c_str() ) == 0 ){  FDRtest = 3; cout << "set BL FDR Method" << endl; }

}

unsigned long int NetworkEnrichment::getSeed(){ return seed; }

bool NetworkEnrichment::_min(double a, double b){ return a<=b ? true : false; }

bool NetworkEnrichment::_max(double a, double b){ return a>=b ? true : false; }
//--------------------------------------------------------


//--------------------------------------------------------
//  UTILITY FUNCTIONS FOR HYPERGEOMETRIC TESTS
//--------------------------------------------------------

//---randomise the distribution of gene id's for each annotation
void NetworkEnrichment::permutation( double tt ){
  
  int k,f,K;

  //--- Initialize random seed:
  // gsl random number generator 'g' in Measures.
  setSeed();

  //linear indexing, size K is rows (N) x cols (F). 
  K=N*F;
  //---reset the permutation matrix
  for(k=0; k<K; k++){ studies[k] = 0; }
  

  for(f=0; f<F; f++){  

    double ANNOTYPESIZE = annoSIZE[f];//annotation size, make sure its been reisze first for the network! -- geneAssociations(), and calculatePvalues()

    int sample_size = 0;

    while( sample_size < ANNOTYPESIZE ){

      int rnd_ind = (int)floor( gsl_rng_uniform(g) * N );

      if( studies[(rnd_ind*F)+f] == 0 ){
	studies[(rnd_ind*F)+f] = 1;
	sample_size++;
      }

    }

  }
    
  
}


//--- Interaction Distance
// See Chapter 5 "Distribution of Interaction Distances",
// Alex T. Kalinka, The probablility of drawing intersections: extending the hypergeometric distribution,
// arXiv:1305.0717v5, (2014).
void NetworkEnrichment::calculateInteractionDistance( int muab, int N, int Cn, int A, int a, int B, int b, int &relD, double &prob ){

  int i,j, D, tempa, tempb;
  vector<pairIntInt> d;
  double mua, mub, prob_A, prob_B;

  relD = 0; prob = 1;
  
  if( a == 0 || b == 0 ) return;

  
  if( b >= a ){ 

    mua = prob_overlap( N, Cn, A, a  );

    for( i=Cn; i>=0; i-- ){
      prob = prob_overlap( N, Cn, B, i );
      if( prob <= mua ){
	tempb=i;
      }
    }

    relD = abs(a-tempb);

  } else {

    mub = prob_overlap( N, Cn, B, b  );

    for( i=Cn; i>=0; i-- ){
      prob = prob_overlap( N, Cn, a, i );
      if( prob <= mub ){
	tempa=i;
      }
    }

    relD = abs(b-tempa);
    
  }
    
  for(i=0; i<=a; i++){
    for(j=0; j<=b; j++){
      //if( abs(i-j) == muab ){ d.push_back( pairIntInt(i,j) ); }
      if( abs(i-j) == relD ){ d.push_back( pairIntInt(i,j) ); }
    }
  }

  prob = 0;
  D    = d.size();  
  for( i=0; i<D; i++ ){
    prob_A = prob_overlap( N, Cn, A, d[i].first  );
    prob_B = prob_overlap( N, Cn, B, d[i].second );    
  
    prob  += prob_A * prob_B;
  }

  if( (prob <= 0) || (prob > 1) ) prob = 1;

}


// The Odds Ratio (OR) and lower / upper 95% confidence interval for 2x2 contingency table
// Szumilas, M. Explaining Odds Ratios, J Can Acad Child Adolesc Psychiatry. 2010 Aug; 19(3): 227–229. 
// Example:
// a = mu
// b = (Cn - mu)
// c = (F  - mu)
// d = (N  - F + mu - Cn)
void NetworkEnrichment::calculateOddsRatio( double a, double b, double c, double d, double &OR, double &lower, double &upper ){

  OR = 0; lower = 0; upper = 0;

  if( (a == 0) || (b == 0) || (c == 0) || (d == 0) ) return;
  
  //the Odds Ratio (OR)
  OR  = (a * d);
  OR /= (b * c);
  //---
  
  //lower 95% confidence interval
  lower = log(OR) - 1.96 * sqrt( (1/a + 1/b + 1/c + 1/d) );

  if( fabs(lower) > MAXEXPO ){
    lower > 0 ? lower = exp( MAXEXPO ) : lower = exp( -MAXEXPO );
  } else {  
    lower = exp( lower ) ;
  }    
  //---
  
  //upper 95% confidence interval
  upper = log(OR) + 1.96 * sqrt( (1/a + 1/b + 1/c + 1/d) );

  if( fabs(upper) > MAXEXPO ){
    upper > 0 ? upper = exp( MAXEXPO ) : upper = exp( -MAXEXPO );
  } else {  
    upper = exp( upper ) ;
  }    
  //---
    
}

//--- FDR Method -> Benjamini and Hochberg FDR (BH), set-up proceedure (1995)
// Benjamini, Y., and Hochberg, Y.
// Controlling the false discovery
// rate:  a practical and powerful approach to multiple testing.
// Journal of the Royal Statistical Society Series B 57
// (1995), 289–300.
void NetworkEnrichment::CalculateFDR_BH( vector<pairDoubInt> order_pv, double _Sigma, double &fdr, double &pv, int &cutoff, vector<pairDoubInt> &_fdr_values ){

  int i,k,N,Kold;
  double _pv, pv_test;
  
  //sort from lowest to highest
  std::sort( order_pv.begin(), order_pv.end(), sortpairDoubInt() );
  
  N = order_pv.size();

  //from lowest to highest
  for(i=0,k=N; i<N; i++, k--){
    pairDoubInt pv_new;
    pv_new   = order_pv[i];
    pv_test  = (double)(_Sigma * (i+1)) / (double) N;
    _pv      = pv_new.first;

    if( _pv <= pv_test ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = k;
    } else {
      break;
    }    

  }
  
  //p-values decreasing, i.e. from max to min p-values
  for(Kold=(N-1),k=(N-1); k>=0; k--){

    pairDoubInt pv_new;
    pairDoubInt pv_old;

    double p_old = 0.0;
    double p_new = 0.0;    
    
    if( (k+1) == N ){
      pv_new        = order_pv[k];
      pv_new.first  = (double)N/(double)(k+1) * pv_new.first;
      _fdr_values.push_back( pairDoubInt(pv_new.first,pv_new.second) );
     } else {
      
      pv_old = order_pv[Kold];
      pv_new = order_pv[k];
      
      p_old  = (double)N/(double)(Kold+1) * pv_old.first;          
      p_new  = (double)N/(double)(k+1)    * pv_new.first;

      if( _min( p_old, p_new ) ){
	_fdr_values.push_back( pairDoubInt(p_old, pv_new.second) );
      } else {
	_fdr_values.push_back( pairDoubInt(p_new, pv_new.second) );
	Kold = k;
      }   

    }
  
    
  }

  for( i=0; i<_fdr_values.size(); i++ ){ _fdr_values[i].first = min( (double)1.0, (double)_fdr_values[i].first ); }
   

}


//--- FDR Method -> Benjamini and Liu FDR (BL), set-down proceedure (1999)
// Benjamini, Y., and Liu, W.
// A step-down multiple hypotheses testing
// procedure  that  controls  the  false  discovery  rate  under  independence.
// Journal of Statistical Planning and Inference 82
// (1999), 163–170
void NetworkEnrichment::CalculateFDR_BL( vector<pairDoubInt> order_pv, double _Sigma, double &fdr, double &pv, int &cutoff, vector<pairDoubInt> &_fdr_values ){

  int i,k,N,Kold;
  double _pv, pv_test;
    
  //sort from lowest to highest
  std::sort( order_pv.begin(), order_pv.end(), sortpairDoubInt() );
  
  N = order_pv.size();

  //from lowest to highest
  for(i=0; i<N; i++){
    pairDoubInt pv_new;
    pv_new   = order_pv[i];
    pv_test  = 1.0 - pow( 1.0 - min(1.0, (double)(_Sigma * N) / (double)( N -(i+1) + 1 )), 1.0/(double)(N-(i+1)+1) );
    _pv      = pv_new.first;

    if( _pv <= pv_test ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = i;
    } else {
      break;
    }    

  }
  
  
  //p-values increasing, i.e. from min to max p-values
  for(Kold=0, k=0; k<N; k++){

    pairDoubInt pv_new;
    pairDoubInt pv_old;

    double p_old = 0.0;
    double p_new = 0.0;

    if( k == 0 ){
      pv_new        = order_pv[k];
      pv_new.first  = ((double)(N-(k+1)+1)/(double)N) * (1.0 - pow((1.0-pv_new.first),(double)(N-(k+1)+1)));
      _fdr_values.push_back( pairDoubInt(pv_new.first, pv_new.second) );
     } else {
      
      pv_old = order_pv[Kold];
      pv_new = order_pv[k];
      
      p_old  = ((double)(N-(Kold+1)+1)/(double)N) * (1.0 - pow((1.0-pv_old.first),(double)(N-(Kold+1)+1)));
      p_new  = ((double)(N-(k+1)+1)/(double)N)    * (1.0 - pow((1.0-pv_new.first),(double)(N-(k+1)+1)));
	
      if( _max( p_old, p_new ) ){
	_fdr_values.push_back( pairDoubInt(p_old, pv_new.second) );
      } else {
	_fdr_values.push_back( pairDoubInt(p_new, pv_new.second) );
	Kold = k;
      }   

    }

    pv_test  = (double)(_Sigma * N) / (double)( N -(k+1) + 1);
    _pv      = pv_new.first;
    
    if( pv_test < _pv ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = (k+1);
    }
    
    
  }
  
  //for( i=0; i<_fdr_values.size(); i++ ){ _fdr_values[i].first = min( (double)1.0, (double)_fdr_values[i].first ); }   

}


//--- FDR Method  -> Benjamini and Yekutieli FDR (BY) set-up proceddure (2001)
//     Benjamini, Y., and Yekutieli, D. (2001).  The control of the false
//     discovery rate in multiple testing under dependency.  _Annals of
//     Statistics_ *29*, 1165-1188.
//     See R's (version 3.4.2) p.adjust function in stats package
void NetworkEnrichment::CalculateFDR_BY( vector<pairDoubInt> order_pv, double _Sigma, double &fdr, double &pv, int &cutoff, vector<pairDoubInt> &_fdr_values ){

  int i,k,N,Kold;
  double q, _pv, pv_test;

  //sort from lowest to highest
  std::sort( order_pv.begin(), order_pv.end(), sortpairDoubInt() );
  
  N = order_pv.size();

  //Normalising constant
  for(q=0.0, i=0; i<N; i++ ){ q += 1.0/(double)(i+1); }

  //from lowest to highest
  for(i=0,k=N; i<N; i++, k--){
    pairDoubInt pv_new;
    pv_new   = order_pv[i];
    pv_test  = (double) (_Sigma * (i+1))/ (double)(q*N);
    _pv      = pv_new.first;
    
    if( _pv <= pv_test ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = k;
    } else {
      break;
    }    

  }
  
  
  //p-values decreasing, i.e. from max to min p-values
  for(Kold=(N-1),k=(N-1); k>=0; k--){

    pairDoubInt pv_new;
    pairDoubInt pv_old;

    double p_old = 0.0;
    double p_new = 0.0;
    
    if( (k+1) == N ){
      pv_new        = order_pv[k];
      pv_new.first  = q * (double)N/(double)(k+1) * pv_new.first;
      _fdr_values.push_back( pairDoubInt(pv_new.first,pv_new.second) );
    } else {
      
      pv_old = order_pv[Kold];
      pv_new = order_pv[k];
      
      p_old  = q * (double)N/(double)(Kold+1) * pv_old.first;          
      p_new  = q * (double)N/(double)(k+1) * pv_new.first;

      if( _min( p_old, p_new ) ){
	_fdr_values.push_back( pairDoubInt(p_old, pv_new.second) );
      } else {
	_fdr_values.push_back( pairDoubInt(p_new, pv_new.second) );
	Kold = k;
      }   

    }    
    
  }

  for( i=0; i<_fdr_values.size(); i++ ){ _fdr_values[i].first = min( (double)1.0, (double)_fdr_values[i].first ); }
  
}



//Hypergeometric distribution.
double NetworkEnrichment::prob_overlap( int Nn, int na, int nb, int nab ){

  double result, temp, num, dem;
  
  result = 0.0; num = 0.0; dem = 0.0; temp = 0.0;

  num = gsl_sf_lnfact( (const unsigned int)na ) + gsl_sf_lnfact( (const unsigned int)(Nn-na) ) + gsl_sf_lnfact( (const unsigned int)nb ) + gsl_sf_lnfact( (const unsigned int)(Nn-nb) );   
  //num = factln( na ) + factln( (Nn-na) )  + factln( nb )  + factln( (Nn-nb) );

  dem = gsl_sf_lnfact( (const unsigned int)Nn ) + gsl_sf_lnfact( (const unsigned int)(na-nab) ) + gsl_sf_lnfact( (const unsigned int)nab ) + gsl_sf_lnfact( (const unsigned int)(Nn-na-nb+nab) ) + gsl_sf_lnfact( (const unsigned int)(nb - nab) );   
  //dem = factln( Nn ) + factln( (na-nab) ) + factln( nab ) + factln( (Nn - na - nb + nab) ) + factln( (nb - nab) );

  temp = num - dem;
  
  if( fabs(temp) > MAXEXPO ){
    temp > 0 ? result = exp( MAXEXPO ) : result = exp( -MAXEXPO );
  } else {  
    result = exp( temp ) ;
  }    
  
  return result;
  
}

// Probability of overlap in network 'N' with three annotation sets: 'A', 'B' and 'C'.
// Alex T. Kalinka, The probablility of drawing intersections: extending the hypergeometric distribution,
// arXiv:1305.0717v5, (2014).
double NetworkEnrichment::prob_overlap( int Nn, int na, int nb, int nc, int v ){

  int i, alpha;
  
  double result, num, dem, temp;
   
  result = 0.0; num = 0.0; dem = 0.0; temp = 0.0;

  alpha = min( (int)(na-v), (int)(nb-v) );  

  for( i=max((int)(na+nb-Nn-v),0); i <= alpha; i++ ){

    if( (int)(Nn-na-nb+v+i) <= 0 || (int)(Nn-nc-i) <= 0 || (int)(Nn-v-i) <= 0 ){
    result += 0.0;
    } else {

      num = gsl_sf_lnfact( (const unsigned int)(na) ) + gsl_sf_lnfact( (const unsigned int)(nb) ) + gsl_sf_lnfact( (const unsigned int)(nc) ) + gsl_sf_lnfact( (const unsigned int)(Nn-na) ) + gsl_sf_lnfact( (const unsigned int)(Nn-nb) ) + gsl_sf_lnfact( (const unsigned int)(Nn-nc) ) +  gsl_sf_lnfact( (const unsigned int)(na-v) ) + gsl_sf_lnfact( (const unsigned int)(Nn-v-i) );
      //num = factln( (int)(na) ) + factln( (int)(nb) ) + factln( (int)(nc) ) + factln( (int)(Nn-na) ) + factln( (int)(Nn-nb) ) + factln( (int)(Nn-nc) ) +  factln( (int)(na-v) ) + factln( (int)(Nn-v-i) );

      dem =  gsl_sf_lnfact( (const unsigned int)(i) ) + gsl_sf_lnfact( (const unsigned int)(v) ) + gsl_sf_lnfact( (const unsigned int)(Nn) ) + gsl_sf_lnfact( (const unsigned int)(Nn) ) +  gsl_sf_lnfact( (const unsigned int)(na-v) ) + gsl_sf_lnfact( (const unsigned int)(nc-v) ) + gsl_sf_lnfact( (const unsigned int)(na-v-i) ) + gsl_sf_lnfact( (const unsigned int)(nb-v-i) ) + gsl_sf_lnfact( (const unsigned int)(Nn-na-nb+v+i) ) + gsl_sf_lnfact( (const unsigned int)(Nn-v-i-nc+v) ); 
      //dem =  factln( (int)(i) ) + factln( (int)(v) ) + factln( (int)(Nn) ) + factln( (int)(Nn) ) +  factln( (int)(na-v) ) + factln( (int)(nc-v) ) + factln( (int)(na-v-i) ) + factln( (int)(nb-v-i) ) + factln( (int)(Nn-na-nb+v+i) ) + factln( (int)(Nn-v-i-nc+v) ); 

      temp = num - dem;
      
      if( fabs(temp) > MAXEXPO ){
	temp > 0 ? result += exp( MAXEXPO ) : result += exp( -MAXEXPO ); 
      } else {  
	result += exp( temp );
      }    

    }

  }

  return result; 
   
}


//is geneID found associated with annotations in annotation set 
void NetworkEnrichment::geneAssociations( int Index, int t, int geneID, int &_assoc, int &_tally ){

  int _t,T;

  T = Fsize[Index];

  //loop over all annotation type 'T'
  for(_t=0; _t<Alines[Index]; _t++){	  
    if( (Alist[Index][_t].ID == geneID) &&	
	(strcmp(Alist[Index][_t].annoID, ANNOS[Index][t].annoID) == 0) ){ 
      _assoc = 1; _tally++; }	    	  
  }
	
}

//is geneID found associated with annotations in annotation set 
void NetworkEnrichment::geneAssociations( int Index, int geneID, int assoc[] ){

  int t,_t,T;

  T = Fsize[Index];

   //loop over all annotation type 'T'
    for(_t=0; _t<Alines[Index]; _t++){	  
      if( (Alist[Index][_t].ID == geneID ) ){	
	
	for(t=0; t<T; t++){	  
	  if( strcmp(Alist[Index][_t].annoID, ANNOS[Index][t].annoID) == 0 ){ 
	    assoc[t] = 1; ANNOS[Index][t].K++;
	  }	    	  
	}
	
      }
    }      

  
}

//Overlap between genes in annotation set and communities
void NetworkEnrichment::overlapinNetwork(){

  int i,m,f,k,K;

  //linear indexing, size K is rows (M) x cols (F).
  K=M*F;

  //--- loop over all communities and each annotation type
  for(k=0; k<K; k++){

    m = floor(k/F);//row index
    f = k % F;     //col index
  
    //--- loop over all genes in the mth cluster which shares the fth annotation type
    for(i=0; i<N; i++){
      if( (geneCOM[i] == (m+1)) && ((int)studies[(i*F)+f] == 1) ){ overlap[(m*F)+f]++; }
    }      
    
  }
  

}

//Overlap between two annotation types in network
void NetworkEnrichment::overlapinNetork( int indexA, int indexB ){

  int i,id,a,b,A,B;

  A = Fsize[indexA];
  B = Fsize[indexB];

  int ind_c [A]; 
  int ind_r [B]; 

  for(a=0; a<A; a++){ ind_c[a]=0; ANNOS[indexA][a].K=0; }//reset
  for(b=0; b<B; b++){ ind_r[b]=0; ANNOS[indexB][b].K=0; }//reset    
    
  //---loop over all ids in network and find overlap between annotation A and annotation type B
  for(i=0; i<N; i++){

    id = Clist[i].ID;//gene ID
      
    for(a=0; a<A; a++){ ind_c[a]=0; }
    for(b=0; b<B; b++){ ind_r[b]=0; }
    
    //loop over all annotation type 'A'
    geneAssociations(indexA, id, ind_c);   

    //loop over all annotation type 'B'
    geneAssociations(indexB, id, ind_r);

    //genes in network and overlap between 'A' and 'B'
    for(a=0; a<A; a++){
      for(b=0; b<B; b++){
	if( ind_c[a] == 1 && ind_r[b] == 1 )
	  overlap[(a*B)+b] += 1;
      }
    }	
    
  }
    
  
}

//Overlap between three annotation types in network
void NetworkEnrichment::overlapinNetork( int indexA, int indexB, int indexC ){

  int i,id,a,b,c,A,B,C;

  A = Fsize[indexA];
  B = Fsize[indexB];
  C = Fsize[indexC];

  int ind_c [A]; 
  int ind_r [B];
  int ind_d [C]; 

  for(a=0; a<A; a++){ ind_c[a]=0; ANNOS[indexA][a].K=0; }//reset
  for(b=0; b<B; b++){ ind_r[b]=0; ANNOS[indexB][b].K=0; }//reset
  for(c=0; c<C; c++){ ind_d[c]=0; ANNOS[indexC][c].K=0; }//reset    
    
  //---loop over all ids in network and find overlap between annotation types A, B and C.
  for(i=0; i<N; i++){

    id = Clist[i].ID;//gene ID
      
    for(a=0; a<A; a++){ ind_c[a]=0; }
    for(b=0; b<B; b++){ ind_r[b]=0; }
    for(c=0; c<C; c++){ ind_d[c]=0; }
    
    //loop over all annotation type 'A'
    geneAssociations(indexA, id, ind_c);   

    //loop over all annotation type 'B'
    geneAssociations(indexB, id, ind_r);

    //loop over all annotation type 'C'
    geneAssociations(indexC, id, ind_d);

    //genes in network and overlap between 'A', 'B' and 'C'
    for(a=0; a<A; a++){
      for(b=0; b<B; b++){
	for(c=0; c<C; c++){
	if( ind_c[a] == 1 && ind_r[b] == 1 && ind_d[c] == 1 )
	  overlap[c+C*(b+B*a)] += 1;
	}
      }	
    }

  }
    
  
}



void NetworkEnrichment::calculateFDR( int Study, int indexA, int indexB, int indexC ){

  int f,m,i,k,a,b,c,A,B,C,K;
  vector<pairDoubInt> pv_sorted;
  vector<pairDoubInt> pv_sortedT;
  vector<pairDoubInt> pv_sortedD;
  vector<pairDoubInt> pv_sortedDT;
  vector<pairDoubInt> pv_sortedRD;

  //reset
  FDR = 0.0; PV = 0.0; LEVEL = 0; TESTS = 0; vector<pairDoubInt> fdr_values;
  vector<pairDoubInt> fdr_valuesD;   vector<pairDoubInt> fdr_valuesT; vector<pairDoubInt> fdr_valuesDT; vector<pairDoubInt> fdr_reldist;
  
  if( Study == 1 ){//calculatePvalues

    //linear indexing, size K is rows (M) x cols (F). 
    K = M*F;
    for(k=0; k<K; k++){
      m = floor(k/F);//row index
      f = k % F;     //col index
      pv_sorted.push_back  ( pairDoubInt(p_values  [(m*F)+f], k) );
      pv_sortedD.push_back ( pairDoubInt(p_valuesD [(m*F)+f], k) );
      pv_sortedDT.push_back( pairDoubInt(p_valuesDT[(m*F)+f], k) );
      pv_sortedT.push_back ( pairDoubInt(p_valuesT [(m*F)+f], k) );
    }    

    //adjust p-values calling 'BY' FDR algorithm 
    if( FDRtest == 1 ){
      CalculateFDR_BY( pv_sorted ,  SIGMA, FDR, PV, LEVEL, fdr_values );
      CalculateFDR_BY( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BY( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BY( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
    }

    //adjust p-values calling 'BH' FDR algorithm 
    if( FDRtest == 2 ){
      CalculateFDR_BH( pv_sorted ,  SIGMA, FDR, PV, LEVEL, fdr_values );
      CalculateFDR_BH( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BH( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BH( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
    }

    //adjust p-values calling 'BL' FDR algorithm 
    if( FDRtest == 3 ){
      CalculateFDR_BL( pv_sorted ,  SIGMA, FDR, PV, LEVEL, fdr_values );
      CalculateFDR_BL( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BL( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BL( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
    }


    //store adjusted pvalues  
    K = fdr_values.size(); 
    TESTS = K;
    for(k=0; k<K; k++){
      int indx  = fdr_values[k].second;
      m = floor(indx/F);//row index
      f = indx % F;     //col index
      padjusted[(m*F)+f]  = fdr_values[k].first;
    }

    K = fdr_valuesD.size(); 
    for(k=0; k<K; k++){    
      int indx = fdr_valuesD[k].second;
      m = floor(indx/F);//row index
      f = indx % F;     //col index      
      padjustedD[(m*F)+f] = fdr_valuesD[k].first;
    }

    K = fdr_valuesDT.size(); 
    for(k=0; k<K; k++){    
      int indx = fdr_valuesDT[k].second;
      m = floor(indx/F);//row index
      f = indx % F;     //col index      
      padjustedDT[(m*F)+f] = fdr_valuesDT[k].first;
    }

    
    K = fdr_valuesT.size(); 
    for(k=0; k<K; k++){    
      int indx = fdr_valuesT[k].second;
      m = floor(indx/F);//row index
      f = indx % F;     //col index      
      padjustedT[(m*F)+f] = fdr_valuesT[k].first;
    }
    
   
  }

  if( Study == 2 ){//calculateOverlapinCommunities

    A = Fsize[indexA];
    B = Fsize[indexB];

    //linear indexing, size K is rows (A*B) x cols (M). 
    K=(A*B)*M;
    for(k=0; k<K; k++){
      i = floor(k/M);//row index
      m = k % M;     //col index
      pv_sorted.push_back  ( pairDoubInt(p_values  [(i*M)+m], k) );
      pv_sortedD.push_back ( pairDoubInt(p_valuesD [(i*M)+m], k) );
      pv_sortedDT.push_back( pairDoubInt(p_valuesDT[(i*M)+m], k) );
      pv_sortedT.push_back ( pairDoubInt(p_valuesT [(i*M)+m], k) );
      if( calRelDist ){
	pv_sortedRD.push_back ( pairDoubInt(p_dist [(i*M)+m], k) );
      }
    }    

    //adjust p-values calling 'BY' FDR algorithm 
    if( FDRtest == 1 ){
      CalculateFDR_BY( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values  );
      CalculateFDR_BY( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BY( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BY( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
      if( calRelDist ){
	CalculateFDR_BY( pv_sortedRD,  SIGMA, FDR, PV, LEVEL, fdr_reldist );
      }
    }

    //adjust p-values calling 'BH' FDR algorithm 
    if( FDRtest == 2 ){
      CalculateFDR_BH( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values   );
      CalculateFDR_BH( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BH( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BH( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
      if( calRelDist ){
	CalculateFDR_BH( pv_sortedRD,  SIGMA, FDR, PV, LEVEL, fdr_reldist );
      }
    }

    //adjust p-values calling 'BL' FDR algorithm 
    if( FDRtest == 3 ){
      CalculateFDR_BL( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values   );
      CalculateFDR_BL( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
      CalculateFDR_BL( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
      CalculateFDR_BL( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
      if( calRelDist ){
	CalculateFDR_BL( pv_sortedRD,  SIGMA, FDR, PV, LEVEL, fdr_reldist );
      }
    }

    
    //store adjusted pvalues  
    K = fdr_values.size(); 
    TESTS = K;
    for(k=0; k<K; k++){
      int indx = fdr_values[k].second;
      i = floor(indx/M);
      m = indx % M;
      padjusted[(i*M)+m] = fdr_values[k].first;
    }

    K = fdr_valuesD.size(); 
    for(k=0; k<K; k++){
      int indx = fdr_valuesD[k].second;
      i = floor(indx/M);
      m = indx % M;
      padjustedD[(i*M)+m] = fdr_valuesD[k].first;
    }

    K = fdr_valuesDT.size(); 
    for(k=0; k<K; k++){
      int indx = fdr_valuesDT[k].second;
      i = floor(indx/M);
      m = indx % M;
      padjustedDT[(i*M)+m] = fdr_valuesDT[k].first;
    }
    
    K = fdr_valuesT.size(); 
    for(k=0; k<K; k++){
      int indx = fdr_valuesT[k].second;
      i = floor(indx/M);
      m = indx % M;
      padjustedT[(i*M)+m] = fdr_valuesT[k].first;
    }

    if( calRelDist ){
      K = fdr_reldist.size(); 
      for(k=0; k<K; k++){
	int indx = fdr_reldist[k].second;
	i = floor(indx/M);
	m = indx % M;
	padjustedRD[(i*M)+m] = fdr_reldist[k].first;
      }

    }
    
  }

  if( Study == 3 ){//calculateOverlapinNetwork

     A = Fsize[indexA];
     B = Fsize[indexB];

     //linear indexing, size K is rows (A) x cols (B). 
     K=A*B;
     for(k=0; k<K; k++){
       a = floor(k/B);//row index
       b = k % B;     //col index
       pv_sorted.push_back  ( pairDoubInt(p_values  [(a*B)+b], k) );
       pv_sortedD.push_back ( pairDoubInt(p_valuesD [(a*B)+b], k) );
       pv_sortedDT.push_back( pairDoubInt(p_valuesDT[(a*B)+b], k) );
       pv_sortedT.push_back ( pairDoubInt(p_valuesT [(a*B)+b], k) );
     }
         
     //adjust p-values calling 'BY' FDR algorithm 
     if( FDRtest == 1 ){
       CalculateFDR_BY( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values );
       CalculateFDR_BY( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BY( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BY( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }

     //adjust p-values calling 'BH' FDR algorithm 
     if( FDRtest == 2 ){
       CalculateFDR_BH( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values );
       CalculateFDR_BH( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BH( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BH( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }
     
     //adjust p-values calling 'BL' FDR algorithm 
     if( FDRtest == 3 ){
       CalculateFDR_BL( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values );
       CalculateFDR_BL( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BL( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BL( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }
     
     //store adjusted pvalues  
     K = fdr_values.size(); 
     TESTS = K;
     for(k=0; k<K; k++ ){
       int indx = fdr_values[k].second;
       a = floor(indx/B);
       b = indx % B;
       padjusted[(a*B)+b] = fdr_values[k].first;
     }

     K = fdr_valuesD.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesD[k].second;
       a = floor(indx/B);
       b = indx % B;
       padjustedD[(a*B)+b] = fdr_valuesD[k].first;
     }

     K = fdr_valuesDT.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesDT[k].second;
       a = floor(indx/B);
       b = indx % B;
       padjustedDT[(a*B)+b] = fdr_valuesDT[k].first;
     }

     
     K = fdr_valuesT.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesT[k].second;
       a = floor(indx/B);
       b = indx % B;
       padjustedT[(a*B)+b] = fdr_valuesT[k].first;
     }

  }

  if( Study == 4 ){//calculateOverlapinNetwork

     A = Fsize[indexA];
     B = Fsize[indexB];
     C = Fsize[indexC];

     //linear indexing, size K is rows (A) x cols (B). 
     K=A*B*C;
     for(k=0; k<K; k++){

       c = k % C;                 //depth index
       b = ((k-c)/C) % B;         //col   index
       a = floor(((k-c)/C - b)/B);//row   index
       
       pv_sorted.push_back  ( pairDoubInt(p_values  [c+C*(b+B*a)], k) );
       pv_sortedD.push_back ( pairDoubInt(p_valuesD [c+C*(b+B*a)], k) );
       pv_sortedDT.push_back( pairDoubInt(p_valuesDT[c+C*(b+B*a)], k) );
       pv_sortedT.push_back ( pairDoubInt(p_valuesT [c+C*(b+B*a)], k) );
     }
         
     //adjust p-values calling 'BY' FDR algorithm 
     if( FDRtest == 1 ){
       CalculateFDR_BY( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values  );
       CalculateFDR_BY( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BY( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BY( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }

     //adjust p-values calling 'BH' FDR algorithm 
     if( FDRtest == 2 ){
       CalculateFDR_BH( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values  );       
       CalculateFDR_BH( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BH( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BH( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }
     
     //adjust p-values calling 'BL' FDR algorithm 
     if( FDRtest == 3 ){
       CalculateFDR_BL( pv_sorted,   SIGMA, FDR, PV, LEVEL, fdr_values  );
       CalculateFDR_BL( pv_sortedD,  SIGMA, FDR, PV, LEVEL, fdr_valuesD );
       CalculateFDR_BL( pv_sortedDT, SIGMA, FDR, PV, LEVEL, fdr_valuesDT );
       CalculateFDR_BL( pv_sortedT,  SIGMA, FDR, PV, LEVEL, fdr_valuesT );
     }
     
     //store adjusted pvalues  
     K = fdr_values.size(); 
     TESTS = K;
     for(k=0; k<K; k++ ){
       int indx = fdr_values[k].second;
       c = indx % C;                 //depth index
       b = ((indx-c)/C) % B;         //col   index
       a = floor(((indx-c)/C - b)/B);//row   index
       padjusted[c+C*(b+B*a)] = fdr_values[k].first;
     }

     K = fdr_valuesD.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesD[k].second;
       c = indx % C;                 //depth index
       b = ((indx-c)/C) % B;         //col   index
       a = floor(((indx-c)/C - b)/B);//row   index
       padjustedD[c+C*(b+B*a)] = fdr_valuesD[k].first;
     }

     K = fdr_valuesDT.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesDT[k].second;
       c = indx % C;                 //depth index
       b = ((indx-c)/C) % B;         //col   index
       a = floor(((indx-c)/C - b)/B);//row   index
       padjustedDT[c+C*(b+B*a)] = fdr_valuesDT[k].first;
     }

     
     K = fdr_valuesT.size(); 
     for(k=0; k<K; k++ ){
       int indx = fdr_valuesT[k].second;
       c = indx % C;                 //depth index
       b = ((indx-c)/C) % B;         //col   index
       a = floor(((indx-c)/C - b)/B);//row   index
       padjustedT[c+C*(b+B*a)] = fdr_valuesT[k].first;
     }
     
     
  }

  
  
  //printFDR();
  
}


//--------------------------------------------------------



//--------------------------------------------------------
//  HYPERGEOMETRIC TEST FUNCTIONS
//--------------------------------------------------------

void NetworkEnrichment::overlapinComsHypergeometricTestRnd( bool runTest ){
  
  int i,m,mm,f,ff,k,K;

  double NORM;

  K    = M*F;
  NORM = K;
   
  
  //--- loop over all communities
  for(m=0; m<M; m++){

    double p_value   = 1.0;
    double p_valueD  = 1.0;
    double p_valueDT = 1.0;
    double p_valueT  = 1.0;
    double Cn        = comSIZE[m];
    
    //continue if community size > MINOVERALP[0]
    if( Cn > MINOVERLAP[0] ){
    
      //--- loop over each annotation type
      for(f=0; f<F; f++ ){

	double study_tally = annoSIZE[f];//K here holds the number of annotation of type f, make sure this is reset for the network rather than the file

	int tally       = 0;
	double mu       = 0;

	//--- loop over all genes in the mth cluster which shares the fth annotation type
	for(i=0; i<N; i++){
	  if( (geneCOM[i] == (m+1)) && ((int)studies[(i*F)+f] == 1) ){ tally++; }
	}

	if( (study_tally > MINOVERLAP[1]) && (tally > MINOVERLAP[2]) ){

	  //calculate p-values
	  p_value   = 0;
	  p_valueD  = 0;
	  p_valueDT = 0;
	  p_valueT  = 0;
	  mu        = 0;
	  mu        = prob_overlap( (int)N, (int)Cn, (int)study_tally, (int)tally );

	  for(i=0; i<=(int)Cn; i++ ){
	    
	    double prob = prob_overlap( (int)N, (int)Cn, (int)study_tally, (int)i );

	    //Enrichment one-sided
	    if( (prob <= mu) && (i <= tally) ){
	      p_value += prob;
	    }

	    //Depletion one-sided
	    if( (prob >= mu) && (i <= tally) ){
	      p_valueD += prob;
	    }
	    
	    //Depletion two-sided
	    if( (prob >= mu) ){
	      p_valueDT += prob;
	    }

	    //Enrichment two-sided
	    if( (prob <= mu) ){
	      p_valueT += prob;
	    }

	  }

	  //---rounding error
	  if( (p_value <= 0)   || (p_value > 1)   ) p_value  = 1.0;

	  if( (p_valueD <= 0)  || (p_valueD > 1)  ) p_valueD = 1.0;

	  if( (p_valueDT <= 0) || (p_valueDT > 1) ) p_valueDT = 1.0;
	   
	  if( (p_valueT <= 0)  || (p_valueT > 1)  ) p_valueT = 1.0;
	     
	} else {
	  p_value   = 1.0;
	  p_valueD  = 1.0;
	  p_valueDT = 1.0;
	  p_valueT  = 1.0;
	}
	
	if( !runTest ){
	  overlap[(m*F)+f]   = (double)tally;

	  //reset
	  p_values[(m*F)+f]   = 0.0;
	  p_valuesD[(m*F)+f]  = 0.0;
	  p_valuesDT[(m*F)+f] = 0.0;
	  p_valuesT[(m*F)+f]  = 0.0;

	  //store p-values from randomised study
	  p_values[(m*F)+f]   = (double)p_value;
	  p_valuesD[(m*F)+f]  = (double)p_valueD;
	  p_valuesDT[(m*F)+f] = (double)p_valueDT;
	  p_valuesT[(m*F)+f]  = (double)p_valueT;
	}

	if( runTest ){
	  //correct for multiple-testing
	  //Test every permute[m][f] value against every pv_values[m][f] value
	  //linear indexing, size K is rows (M) x cols (F).
	  for(k=0; k<K; k++){
	    mm = floor(k/F);
	    ff = k % F;
	    if( (double)p_value <= fabs(p_values[(mm*F)+ff]) )
	      { permute[(mm*F)+ff] = permute[(mm*F)+ff] + 1/NORM; }

	    if( (double)p_valueD <= fabs(p_valuesD[(mm*F)+ff]) )
	      { permuteD[(mm*F)+ff] = permuteD[(mm*F)+ff] + 1/NORM; }

	    if( (double)p_valueDT <= fabs(p_valuesDT[(mm*F)+ff]) )
	      { permuteDT[(mm*F)+ff] = permuteDT[(mm*F)+ff] + 1/NORM; }
	    
	    if( (double)p_valueT <= fabs(p_valuesT[(mm*F)+ff]) )
	      { permuteT[(mm*F)+ff] = permuteT[(mm*F)+ff] + 1/NORM; }	    
	  }
	}//runTest
      }//F
      
    } else {

      if( runTest ){
	//--- loop over each annotation type
	for(f=0; f<F; f++ ){		
	  //correct for multiple-testing
	  //Test every permute[m][f] value against every pv_values[m][f] value
	  //linear indexing, size K is rows (M) x cols (F).
	  for(k=0; k<K; k++){
	    mm = floor(k/F);
	    ff = k % F;
	    if( (double)p_value <= fabs(p_values[(mm*F)+ff]) )
	      { permute[(mm*F)+ff] = permute[(mm*F)+ff] + 1/NORM; }

	    if( (double)p_valueD <= fabs(p_valuesD[(mm*F)+ff]) )
	      { permuteD[(mm*F)+ff] = permuteD[(mm*F)+ff] + 1/NORM; }

	    if( (double)p_valueDT <= fabs(p_valuesDT[(mm*F)+ff]) )
	      { permuteDT[(mm*F)+ff] = permuteDT[(mm*F)+ff] + 1/NORM; }
	    
	    if( (double)p_valueT <= fabs(p_valuesT[(mm*F)+ff]) )
	      { permuteT[(mm*F)+ff] = permuteT[(mm*F)+ff] + 1/NORM; }	    

	  }	
	}
      }//runTest
    }//if	  
    
  }//M
 
}


void NetworkEnrichment::overlapinComsHypergeometricTest(){

  
  int i,m,f,k,K;

  //linear indexing, size K is rows (M) x cols (F).
  K = M*F;
  //---loop over all communities (M) and each annotation type (F) 
  for(k=0; k<K; k++){

    m = floor(k/F);//row index
    f = k % F;     //col index

    double Cn = comSIZE[m];
    double study_tally = annoSIZE[f];//ANNOS[f].K;//K here holds the number of annotation of type f, make this is reset for the network rather than the file

    double p_value   = 0.0;
    double p_valueD  = 0.0;
    double p_valueDT = 0.0;
    double p_valueT  = 0.0;
    double mu        = 0;
    double tally     = overlap[(m*F)+f];

    //continue if no: annotation of type f > MINOVERALP[1] and
    //         if overlap of annotation of type f in community m > MINOVERLAP[2]
    if( (Cn > MINOVERLAP[0] ) && (study_tally > MINOVERLAP[1]) && (tally > MINOVERLAP[2]) ){

      //calculate p-values
      p_value   = 0;
      p_valueD  = 0;
      p_valueDT = 0;
      p_valueT  = 0;
      mu        = 0;
      mu        = prob_overlap( (int)N, (int)Cn, (int)study_tally, (int)tally );

      for(i=0; i<=(int)Cn; i++ ){

	double prob = prob_overlap( (int)N, (int)Cn, (int)study_tally, (int)i );

	//Enrichment one-sided
	if( (prob <= mu) && (i <= tally) ){
	  p_value += prob;
	}

	//Depletion one-sided
	if( (prob >= mu) && (i <= tally) ){
	  p_valueD += prob;
	}

	//Depletion two-sided
	if( (prob >= mu) ){
	  p_valueDT += prob;
	}

	//Enrichment two-sided
	if( (prob <= mu) ){
	  p_valueT += prob;
	}


      }

      //---rounding error
      if( (p_value <= 0)   || (p_value > 1)  ) p_value  = 1.0;

      if( (p_valueD <= 0)  || (p_valueD > 1) ) p_valueD = 1.0;

      if( (p_valueDT <= 0) || (p_valueDT > 1)) p_valueDT = 1.0;

      if( (p_valueT <= 0)  || (p_valueT > 1) ) p_valueT = 1.0;
      

      //calculate p-values, enrichment one-side
      p_values[(m*F)+f]   = (double)p_value;

      //calculate p-values, depletion
      p_valuesD[(m*F)+f]  = (double)p_valueD;

      //calculate p-values, depletion
      p_valuesDT[(m*F)+f] = (double)p_valueDT;

      //calculate p-values, enrichment two-sided
      p_valuesT[(m*F)+f]  = (double)p_valueT;
      
    } else {

      p_values[(m*F)+f]   = 1.0;

      p_valuesD[(m*F)+f]  = 1.0;

      p_valuesDT[(m*F)+f] = 1.0;

      p_valuesT[(m*F)+f]  = 1.0;
      
    }
        
  }
   
  
}

void NetworkEnrichment::overlapinComsHypergeometricTest(int indexA, int indexB){

  int i,k,m,a,b,As,Bs;

  As = Fsize[indexA];
  Bs = Fsize[indexB];  
  
  //--- loop over each Annotation type
  for(k=0, a=0; a<As; a++){
    for(b=0; b<Bs; b++, k++){
       
      double tot = 0.0;
      double A   = ANNOS[indexA][a].K;
      double B   = ANNOS[indexB][b].K;
      
      //--- loop over all communities
      for(m=0; m<M; m++){

	double p_value   = 1.0;
	double p_valueD  = 1.0;
	double p_valueDT = 1.0;
	double p_valueT  = 1.0;
	double mu        = 0.0;
	double prob_A    = 0.0;
	double prob_B    = 0.0;

	int tally        = 0;
	int tally_na     = 0;
	int tally_nb     = 0;

	
	//--- loop over all genes in the mth cluster which shares the fth Disease type
	for(i=0; i<N; i++){
	  
	  if( Clist[i].K == (m+1) ){
	    
	    int temp_na = -1;
	    int temp_nb = -1;
	    
	    //loop over all annotation type 'A'
	    geneAssociations( indexA, a, Clist[i].ID, temp_na, tally_na );
	    
	    //loop over all annotation type 'B'
	    geneAssociations( indexB, b, Clist[i].ID, temp_nb, tally_nb );
	    
	    if( temp_na == 1 && temp_nb == 1 ){ tally++; }
	    
	  }
	}//network ids	 	
	

	//overlap of annotation types A and B in community m
	muCab[(k*M)+m] = tally;

	//product of annotation types A and B in community m
	nab[(k*M)+m] = tally_na * tally_nb;
	
	if( (comSIZE[m] > MINOVERLAP[0] ) && (overlap[(a*Bs)+b] > MINOVERLAP[1]) && (tally > MINOVERLAP[2]) ){ 

	  p_value   = 0;
	  p_valueD  = 0;
	  p_valueDT = 0;
	  p_valueT  = 0;
	  mu        = 0;
	  prob_A    = prob_overlap( (int)N, (int)comSIZE[m], (int)A, (int)tally_na );
	  prob_B    = prob_overlap( (int)N, (int)comSIZE[m], (int)B, (int)tally_nb );

	  mu        = prob_overlap( (int)N, (int)comSIZE[m], (int)overlap[(a*Bs)+b], (int)tally );
	  mu       *= (prob_A * prob_B);

	  if( calRelDist ){
	    double prob_RD = 0.0;
	    int    relDist = 0;
	    
	    calculateInteractionDistance( (int)overlap[(a*Bs)+b], (int)N, (int)comSIZE[m], (int)A, (int)tally_na, (int)B, (int)tally_nb, relDist, prob_RD );

	    p_dist[(k*M)+m]  = prob_RD;
	    reldist[(k*M)+m] = relDist;
	  }
	    
	   for(i=0; i<=(int)comSIZE[m]; i++ ){

	     double prob  = prob_overlap( (int)N, (int)comSIZE[m], (int)overlap[(a*Bs)+b], (int)i );
	            prob *= (prob_A * prob_B);

	     //Enrichment one-sided
	     if( (prob <= mu) && (i <= tally) ){
	       p_value += prob;
	     }

	     //Depletion one-sided
	     if( (prob >= mu) && (i <= tally) ){
	       p_valueD += prob;
	     }

	     //Depletion two-sided
	     if( (prob >= mu) ){
	       p_valueDT += prob;
	     }

	     //Enrichment, two-sided
	     if( (prob <= mu) ){
	       p_valueT += prob;
	     }

	   }
	     
	   //---rounding error
	   if( (p_value <= 0.0)   || (p_value > 1.0)   ) p_value  = 1.0;

	   if( (p_valueD <= 0.0)  || (p_valueD > 1.0)  ) p_valueD = 1.0;

	   if( (p_valueDT <= 0.0) || (p_valueDT > 1.0) ) p_valueDT = 1.0;

	   if( (p_valueT <= 0.0)  || (p_valueT > 1.0)  ) p_valueT = 1.0;
	     
	} else {
	     
	       p_values[(k*M)+m]   = 1.0;

	       p_valuesD[(k*M)+m]  = 1.0;

	       p_valuesDT[(k*M)+m] = 1.0;

	       p_valuesT[(k*M)+m]  = 1.0;

	}

	//calculate p-values, enrichment one-sided
	p_values[(k*M)+m]   = (double)p_value;

	//calculate p-values, depletion
	p_valuesD[(k*M)+m]  = (double)p_valueD;

	//calculate p-values, depletion
	p_valuesDT[(k*M)+m] = (double)p_valueDT;
	
	//calculate p-values, enrichment two-sided
	p_valuesT[(k*M)+m]  = (double)p_valueT;
	
      }
    }
  }	     

}

void NetworkEnrichment::overlapinNetHypergeometricTest(int indexA, int indexB){

  int i,k,K,a,b,A,B;

  //indexA ==> ROW index
  //indexB ==> COL index
  A = Fsize[indexA];
  B = Fsize[indexB];

  //calculate muab
  //linear indexing, size K is rows (A) x cols (B).
  K=A*B;
  for(k=0; k<K; k++){
    
    a = floor(k/B);//row index
    b = k % B;     //col index

    muab[(a*B)+b] = 0.0;
    muab[(a*B)+b] = prob_overlap( (int)N, (int)ANNOS[indexA][a].K, (int)ANNOS[indexB][b].K, (int)overlap[(a*B)+b] );  
  } 

  //---calculate prob
  //linear indexing, size K is rows (A) x cols (B).
  for(k=0; k<K; k++){
    
    a = floor(k/B);
    b = k % B;
    
    double p_value   = 0.0;
    double p_valueD  = 0.0;
    double p_valueDT = 0.0;
    double p_valueT  = 0.0;
    double MIN       = 0.0;

    //---if overlap <= MINOVERLAP, overlap too small
    if( (double)overlap[(a*B)+b]   > MINOVERLAP[0] &&
	(double)ANNOS[indexA][a].K > MINOVERLAP[1] &&
	(double)ANNOS[indexB][b].K > MINOVERLAP[2] ){

      MIN = (double)ANNOS[indexA][a].K;
      if( (double)ANNOS[indexB][b].K < MIN ){ MIN = (double)ANNOS[indexB][b].K; }

      for(i=0; i<=MIN; i++){
	double prob = prob_overlap( (int)N, (int)ANNOS[indexA][a].K, (int)ANNOS[indexB][b].K, (int)i );

	//Enrichment one-sided
	if( (prob <= muab[(a*B)+b]) && (i <= overlap[(a*B)+b]) ){
	  p_value += prob;
	}

	//Depletion one-sided
	if( (prob > muab[(a*B)+b]) && (i <= overlap[(a*B)+b]) ){
	  p_valueD += prob;
	}

	//Depletion two-sided
	if( (prob > muab[(a*B)+b]) ){
	  p_valueDT += prob;
	}

	//Enrichment two-sided
	if( (prob <= muab[(a*B)+b]) ){
	  p_valueT += prob;
	}

      }

    } else {
      p_value   = 1.0;
      p_valueD  = 1.0;
      p_valueDT = 1.0;
      p_valueT  = 1.0;
    }

     //---rounding error
    if( (p_value <= 0)   || (p_value > 1)   ) p_value  = 1.0;

    if( (p_valueD <= 0)  || (p_valueD > 1)  ) p_valueD = 1.0;

    if( (p_valueDT <= 0) || (p_valueDT > 1) ) p_valueDT = 1.0;

    if( (p_valueT <= 0)  || (p_valueT > 1)  ) p_valueT = 1.0;
    

    p_values[(a*B)+b]   = p_value;

    p_valuesD[(a*B)+b]  = p_valueD;

    p_valuesDT[(a*B)+b] = p_valueDT;

    p_valuesT[(a*B)+b]  = p_valueT;

  }

 
}

void NetworkEnrichment::overlapinNetHypergeometricTest(int indexA, int indexB, int indexC){

  int i,k,K,a,b,c,A,B,C;

  //indexA ==> ROW index
  //indexB ==> COL index
  //indexC ==> DEP index
  A = Fsize[indexA];
  B = Fsize[indexB];
  C = Fsize[indexC];

  //calculate muab
  //linear indexing, size K is rows (A) x cols (B) x depth (C).
  K=A*B*C;
  for(k=0; k<K; k++){

    c = k % C;                 //depth index
    b = ((k-c)/C) % B;         //col   index
    a = floor(((k-c)/C - b)/B);//row   index
    
    muab[c+C*(b+B*a)] = 0.0;
    muab[c+C*(b+B*a)] = prob_overlap( (int)N, (int)ANNOS[indexA][a].K, (int)ANNOS[indexB][b].K, (int)ANNOS[indexC][c].K, (int)overlap[c+C*(b+B*a)] );  
      
  } 

  //---calculate prob
  //linear indexing, size K is rows (A) x cols (B) x depth (C).
  for(k=0; k<K; k++){
    
    c = k % C;                 //depth index
    b = ((k-c)/C) % B;         //col   index
    a = floor(((k-c)/C - b)/B);//row   index
    
    double p_value   = 0.0;
    double p_valueD  = 0.0;
    double p_valueDT = 0.0;
    double p_valueT  = 0.0;
    double MIN       = 0.0;
    
    //---if overlap <= MINOVERLAP, overlap too small
    if( (double)overlap[c+C*(b+B*a)] > MINOVERLAP[0] &&
	(double)ANNOS[indexA][a].K   > MINOVERLAP[1] &&
	(double)ANNOS[indexB][b].K   > MINOVERLAP[2] &&
	(double)ANNOS[indexC][c].K   > MINOVERLAP[3] ){

      if( (ANNOS[indexA][a].K <= ANNOS[indexB][b].K) && (ANNOS[indexA][a].K <= ANNOS[indexC][c].K) )
	MIN = ANNOS[indexA][a].K;
      else if (  (ANNOS[indexB][b].K <= ANNOS[indexA][a].K) && (ANNOS[indexB][b].K <= ANNOS[indexC][c].K) )
	MIN = ANNOS[indexB][b].K;
      else
	MIN = ANNOS[indexC][c].K;

      for(i=0; i<=MIN; i++){
	
	double prob = prob_overlap( (int)N, (int)ANNOS[indexA][a].K, (int)ANNOS[indexB][b].K, (int)ANNOS[indexC][c].K, (int)i );

	//Enrichment one-sided
	if( (prob <= muab[c+C*(b+B*a)]) && (i <= overlap[c+C*(b+B*a)]) ){
	  p_value += prob;
	}

	//Depletion one-sided
	if( (prob > muab[c+C*(b+B*a)]) && (i <= overlap[c+C*(b+B*a)]) ){
	  p_valueD += prob;
	}

	//Depletion two-sided
	if( (prob > muab[c+C*(b+B*a)]) ){
	  p_valueDT += prob;
	}

	//Enrichmnet two-sided
	if( (prob <= muab[c+C*(b+B*a)]) ){
	  p_valueT += prob;
	}

      }

    } else {
      p_value   = 1.0;
      p_valueD  = 1.0;
      p_valueDT = 1.0;
      p_valueT  = 1.0;
    }

    //---rounding error
    if( (p_value <= 0)   || (p_value > 1)   ) p_value  = 1.0;

    if( (p_valueD <= 0)  || (p_valueD > 1)  ) p_valueD = 1.0;

    if( (p_valueDT <= 0) || (p_valueDT > 1) ) p_valueDT = 1.0;
    
    if( (p_valueT <= 0)  || (p_valueT > 1)  ) p_valueT = 1.0;
    
    p_values[c+C*(b+B*a)]   = p_value;

    p_valuesD[c+C*(b+B*a)]  = p_valueD;

    p_valuesDT[c+C*(b+B*a)] = p_valueDT;

    p_valuesT[c+C*(b+B*a)]  = p_valueT;


  }
  
  
}
//--------------------------------------------------------



//--------------------------------------------------------
//  PRINT FUNCTIONS
//--------------------------------------------------------
void NetworkEnrichment::printFDR(){

  cout << "" << endl;
  cout << "-----------------------" << endl;
  cout << "FDR method: " << FDRmethods[FDRtest-1] << endl;
  cout << "Sig. level: " << SIGMA << endl;
  cout << "No: TESTS : " << TESTS << endl;
  cout << "LEVEL     : " << LEVEL << endl;
  cout << "-----------------------" << endl;
  cout << " p.values <= " << PV << " sig. at test statistic: " << FDR << endl;
  cout << "-----------------------" << endl;
  
}

void NetworkEnrichment::printOverlapinCommunities( const char *outdir, const char *ext, bool printFDR, bool printPerm ){

  
  int i,m,f,k,K;
  
  for(i=0; i<SIGMASIZE; i++){
    bonferroni[i] = sigma[i]/(double)(M*F);//Bonferroni at sigma[i]
  }

  
  char buffer    [BUFFERSIZE];
  char bufferOR  [BUFFERSIZE];
  
  sprintf(buffer ,"%s/%s_%s.csv",outdir,baseNAME[0].c_str(),ext);
  fstream *fileout = new fstream(buffer,ios_base::out);

  (*fileout)  << "C" << dels[0] << "Cn" << dels[0];

    for(f=0; f<F; f++){
      if( printALT ){
	(*fileout)  << (printID ? ANNOS[ANNOindex][f].annoID : ANNOS[ANNOindex][f].annoDES) << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "OR" << dels[0] << "95% CI" << dels[0] << "p-value" << dels[0] << "p.adjusted" << dels[0] << "{p_value}" << dels[0] << "BC" << dels[0] << "p-value(ALT)" << dels[0] << "p.adjusted(ALT)" << dels[0] << "BC(ALT)" << dels[0];
      } else {
	(*fileout)  << (printID ? ANNOS[ANNOindex][f].annoID : ANNOS[ANNOindex][f].annoDES) << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "OR" << dels[0] << "95% CI" << dels[0] << "p-value" << dels[0] << "p.adjusted" << dels[0] << "{p_value}" << dels[0] << "BC" << dels[0];
      }

    }

    //--- loop over all communities
    for(m=0; m<M; m++){

      (*fileout) << "" << endl;
      (*fileout) << (COMS[m].first-KOFFSET) << dels[0] << COMS[m].second << dels[0];

      //--- loop over each annotation type
      for(f=0; f<F; f++ ){
      
	string starsEO = "";
	string starsET = "";	
	string starsDO = "";
	string starsDT = "";
	
	for(i=0; i<SIGMASIZE; i++){
	  if( p_values[(m*F)+f]   <= bonferroni[i] ){ starsEO += "*";  }
	  if( p_valuesT[(m*F)+f]  <= bonferroni[i] ){ starsET += "*";  }
	  if( p_valuesD[(m*F)+f]  <= bonferroni[i] ){ starsDO += "*";  }
	  if( p_valuesDT[(m*F)+f] <= bonferroni[i] ){ starsDT += "*";  }
	}

	//--- fold-change
	//double fc = 0.0;
	//fc        = double(overlap[(m*F)+f]) / double(COMS[m].second);
	//fc       /= double(COMS[m].second) / double(N);

	//--- Hypergeometric mean
	double mn = 0.0;
	mn        = double(COMS[m].second) * double(ANNOS[ANNOindex][f].K);
	mn       /= double(N);


	//--- Hypergeometric variance
	//double vr = 0.0;
	//vr        = mn * double(N-ANNOS[ANNOindex][f].K) * double(N-COMS[m].second);
	//vr       /= double(N) * double(N-1);
	
	//--- Odds Ratio
	double ors = 0;
	double orL = 0;
	double orU = 0;

	double a   = double(overlap[(m*F)+f]);
	double b   = double(COMS[m].second - overlap[(m*F)+f]);
	double c   = double(ANNOS[ANNOindex][f].K - overlap[(m*F)+f]);
	double d   = double(N - ANNOS[ANNOindex][f].K + overlap[(m*F)+f] - COMS[m].second);
	
	calculateOddsRatio( a, b, c, d, ors, orL, orU );

	sprintf(bufferOR,"[%.2f,%.2f]",orL, orU);	
	//---
	
	//--- print Anno size
	sprintf(buffer,"%d",ANNOS[ANNOindex][f].K);
	
	if( printALT ){
	  (*fileout) << (printAn ? string(buffer) : "") << dels[0] << overlap[(m*F)+f] << dels[0] << mn << dels[0] << ors << dels[0] << bufferOR << dels[0] << (printTwoSided ? p_valuesT[(m*F)+f] : p_values[(m*F)+f]) << dels[0] << (printTwoSided ? padjustedT[(m*F)+f] : padjusted[(m*F)+f]) << dels[0] << ( (this->pesudocount+permute[(m*F)+f])/NoP) * 100 << dels[0] << (printTwoSided ? starsET : starsEO) << dels[0] << (printTwoSided ? p_valuesDT[(m*F)+f] : p_valuesD[(m*F)+f]) << dels[0] << (printTwoSided ? padjustedDT[(m*F)+f] : padjustedD[(m*F)+f]) << dels[0] << (printTwoSided ? starsDT : starsDO) << dels[0];
	} else {
	  (*fileout) << (printAn ? string(buffer) : "") << dels[0] << overlap[(m*F)+f] << dels[0] << mn << dels[0] << ors << dels[0] << bufferOR << dels[0] << (printTwoSided ? p_valuesT[(m*F)+f] : p_values[(m*F)+f]) << dels[0] << (printTwoSided ? padjustedT[(m*F)+f] : padjusted[(m*F)+f]) << dels[0] << ( (this->pesudocount+permute[(m*F)+f])/NoP) * 100 << dels[0] << (printTwoSided ? starsET : starsEO) << dels[0];
	}	
      }
    }
      
    (*fileout) << "" << endl;
    
    fileout->close();
   
   
}

void NetworkEnrichment::printOverlapinCommunitiesAlt( const char *outdir, const char *ext, bool printFDR){

  int i,m,f,k,K;

  for(i=0; i<SIGMASIZE; i++){
    bonferroni[i] = sigma[i]/(double)(M*F);//Bonferroni at sigma[i]
  }

  
  char buffer    [BUFFERSIZE];
  char bufferOR  [BUFFERSIZE];

  sprintf(buffer ,"%s/%s_%s.csv",outdir,baseNAME[3].c_str(),ext);
  fstream *fileout = new fstream(buffer,ios_base::out);

  (*fileout) << "N: " << N << endl;
  if( printALT ){
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "OR" << dels[0] << "95% CI" << dels[0] << "p-value" << dels[0] << "p.adjusted" << dels[0] << "BC" << dels[0] << "p-value(ALT)" << dels[0] << "p.adjusted(ALT)" << dels[0] << "BC" << endl;
  } else {  
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0]  << "OR" << dels[0] << "95% CI" << "p-value" << dels[0] << "p.adjusted" << dels[0] << "BC" << endl;
  }

  //linear indexing, size K is rows (M) x cols (F). 
  K=M*F;
  for(k=0; k<K; k++){

    m = floor(k/F);//row index
    f = k % F;     //col index

    string starsEO = "";
    string starsET = "";	
    string starsDO = "";
    string starsDT = "";
	
    for(i=0; i<SIGMASIZE; i++){
      if( p_values[(m*F)+f]   <= bonferroni[i] ){ starsEO += "*";  }
      if( p_valuesT[(m*F)+f]  <= bonferroni[i] ){ starsET += "*";  }
      if( p_valuesD[(m*F)+f]  <= bonferroni[i] ){ starsDO += "*";  }
      if( p_valuesDT[(m*F)+f] <= bonferroni[i] ){ starsDT += "*";  }
    }

    sprintf(buffer,"community_%d",(COMS[m].first-KOFFSET));

    //--- fold-chnage
    //double fc = 0.0;
    //fc        = double(overlap[(m*F)+f]) / double(COMS[m].second);
    //fc       /= double(COMS[m].second) / double(N);

    //--- Hypergeometric mean, or expected overlap
    double mn = 0.0;
    mn        = double(COMS[m].second) * double(ANNOS[ANNOindex][f].K);
    mn       /= double(N);

    //--- Hypergeometric variance
    //double vr = 0.0;
    //vr        = mn * double(N-ANNOS[ANNOindex][f].K) * double(N-COMS[m].second);
    //vr       /= double(N) * double(N-1);

    //--- Odds Ratio
    double ors = 0;
    double orL = 0;
    double orU = 0;

    double a   = double(overlap[(m*F)+f]);
    double b   = double(COMS[m].second - overlap[(m*F)+f]);
    double c   = double(ANNOS[ANNOindex][f].K - overlap[(m*F)+f]);
    double d   = double(N - ANNOS[ANNOindex][f].K + overlap[(m*F)+f] - COMS[m].second);
	
    calculateOddsRatio( a, b, c, d, ors, orL, orU );

    sprintf(bufferOR,"[%.2f,%.2f]",orL, orU);	
    //---
	
    
    
    if( printALT ){
      (*fileout) << string(buffer) << dels[0] << COMS[m].second << dels[0] << (printID ? ANNOS[ANNOindex][f].annoID : ANNOS[ANNOindex][f].annoDES) << dels[0] << ANNOS[ANNOindex][f].K << dels[0] << overlap[(m*F)+f] << dels[0] << mn << dels[0] << ors << dels[0] << bufferOR << dels[0] << (printTwoSided ? p_valuesT[(m*F)+f] : p_values[(m*F)+f]) << dels[0] << (printTwoSided ?  padjustedT[(m*F)+f] : padjusted[(m*F)+f]) << dels[0] << (printTwoSided ? starsET : starsEO) << dels[0] << (printTwoSided ? p_valuesDT[(m*F)+f] : p_valuesD[(m*F)+f]) << dels[0] << (printTwoSided ? padjustedDT[(m*F)+f] : padjustedD[(m*F)+f]) << dels[0] << (printTwoSided ? starsDT : starsDO) << endl;
    } else {      
      (*fileout) << string(buffer) << dels[0] << COMS[m].second << dels[0] << (printID ? ANNOS[ANNOindex][f].annoID : ANNOS[ANNOindex][f].annoDES) << dels[0] << ANNOS[ANNOindex][f].K << dels[0] << overlap[(m*F)+f] << dels[0] << mn << dels[0] << ors << dels[0] << bufferOR << dels[0] << (printTwoSided ? p_valuesT[(m*F)+f] : p_values[(m*F)+f]) << dels[0] << (printTwoSided ?  padjustedT[(m*F)+f] : padjusted[(m*F)+f]) << dels[0] << (printTwoSided ? starsET : starsEO) << endl;
      }

  } 

  (*fileout) << "" << endl;
  
  fileout->close();
 
   
}


void NetworkEnrichment::printOverlapinCommunities( int indexA, int indexB, const char *outdir, const char *ext, bool printFDR){
  
  int i,k,m,a,b,A,B;

  A = Fsize[indexA];
  B = Fsize[indexB]; 
  
  for(i=0; i<SIGMASIZE; i++){
    bonferroni[i] = sigma[i]/(double)(A*B*M);//Bonferroni at sigma[i]
  }
  
  char buffer   [BUFFERSIZE];
  char buffer1  [BUFFERSIZE];
  char buffer2  [BUFFERSIZE];
  char buffer3  [BUFFERSIZE];
  char buffer4  [BUFFERSIZE];
  sprintf(buffer ,"%s/%s_%s.csv",outdir,baseNAME[1].c_str(),ext);
  fstream *fileout = new fstream(buffer,ios_base::out);

  (*fileout) << "" << dels[0] << "" << dels[0] << "" << dels[0] << "" << dels[0] << "Communities" << dels[0];
  
  //--- loop over all communities
  for(m=0; m<M; m++){
    if( printALT ){ (*fileout) << (m+1-KOFFSET) << dels[0] << "" << dels[0] << "" << dels[0]; }
    else { (*fileout) << (m+1-KOFFSET) << dels[0] << "" << dels[0]; }
  }

  (*fileout) << "" << endl;    
  (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "actual overlap" << dels[0];
  
  //--- loop over all communities
  for(m=0; m<M; m++){
    if( printALT ){ (*fileout) << "{a.overlap,e.overlap,OR,95% CI}" << dels[0] << "{p.value,p.adjusted,BC}" << dels[0] << "{p.value(ALT),p.adjusted(ALT),BC(ALT)}" << dels[0]; } else {
      (*fileout) << "{a.overlap,e.overlap,OR,95% CI}" << dels[0] << "{p.value,p.adjusted,BC}" << dels[0]; }

    if( calRelDist ){
      (*fileout) << "{relDist, p.value, p.adjusted}" << dels[0];
    }
  }

  (*fileout) << "" << endl;

  //--- loop over each Annotation type
  for(k=0, a=0; a<A; a++){
    for(b=0; b<B; b++, k++){
      (*fileout) <<  (printID ? ANNOS[indexA][a].annoID : ANNOS[indexA][a].annoDES) << dels[0] << ANNOS[indexA][a].K << dels[0] <<  (printID ? ANNOS[indexB][b].annoID : ANNOS[indexB][b].annoDES) << dels[0] << ANNOS[indexB][b].K << dels[0] << overlap[(a*B)+b] << dels[0];

      //--- loop over all communities
      for(m=0; m<M; m++){

	string starsEO = "";
	string starsET = "";	
	string starsDO = "";
	string starsDT = "";
	
	for(i=0; i<SIGMASIZE; i++){
	  if( p_values  [(k*M)+m] <= bonferroni[i] ){ starsEO += "*";  }
	  if( p_valuesT [(k*M)+m] <= bonferroni[i] ){ starsET += "*";  }
	  if( p_valuesD [(k*M)+m] <= bonferroni[i] ){ starsDO += "*";  }
	  if( p_valuesDT[(k*M)+m] <= bonferroni[i] ){ starsDT += "*";  }
	}
	
	//--- fold-chnage
	//double fc = 0.0;
	//fc        = double(muCab[(k*M)+m]) / double(COMS[m].second);
	//fc       /= double(COMS[m].second) / double(N);

	//--- mean
	double mn = 0.0;
	mn        = double(COMS[m].second) * double(overlap[(a*B)+b]);
	mn       /= double(N);

	//--- Odds Ratio
	double ors = 0;
	double orL = 0;
	double orU = 0;

	double ORa = double(muCab[(k*M)+m]);
	double ORb = double(COMS[m].second - muCab[(k*M)+m]);
	double ORc = double(overlap[(a*B)+b] - muCab[(k*M)+m]);
	double ORd = double(N - overlap[(a*B)+b] + muCab[(k*M)+m] - COMS[m].second);
	
	calculateOddsRatio( ORa, ORb, ORc, ORd, ors, orL, orU );
	//---

	if( calRelDist ){
	  sprintf(buffer4,"%.1f, %G, %G", reldist[(k*M)+m], p_dist[(k*M)+m], padjustedRD[(k*M)+m]);
	}
	  
	sprintf(buffer1,"%.1f, %.2f, %.2f, [%.2f,%.2f], %G",muCab[(k*M)+m],mn,ors,orL, orU);

	if( printTwoSided ){
	  sprintf(buffer2,"%G, %G, %s",p_valuesT [(k*M)+m],padjustedT [(k*M)+m],starsET.c_str());
	  sprintf(buffer3,"%G, %G, %s",p_valuesDT[(k*M)+m],padjustedDT[(k*M)+m],starsDT.c_str());
	} else {
	  sprintf(buffer2,"%G, %G, %s",p_values [(k*M)+m],padjusted [(k*M)+m],starsEO.c_str());
	  sprintf(buffer3,"%G, %G, %s",p_valuesD[(k*M)+m],padjustedD[(k*M)+m],starsDO.c_str());
	}

	if( printALT ){
	  (*fileout) << buffer1 << dels[0] << buffer2 << dels[0] << buffer3 << dels[0]; 
	} else {
	  (*fileout) << buffer1 << dels[0] << buffer2 << dels[0];
	}

	if( calRelDist ){
	  (*fileout) << buffer4 << dels[0];
	}
	
      }      
      (*fileout) << "" << endl;
    }
  }

  (*fileout) << "" << endl;
  
  fileout->close();
      
}

void NetworkEnrichment::printOverlapinNetwork( int indexA, int indexB, const char *outdir, const char *ext, bool printFDR){
  
  int i,k,K,a,b,A,B;

  A = Fsize[indexA];
  B = Fsize[indexB]; 
  
  for(i=0; i<SIGMASIZE; i++){
    bonferroni[i] = sigma[i]/(double)(A*B);//Bonferroni at sigma[i]
  }

  
  char buffer  [BUFFERSIZE];
  sprintf(buffer ,"%s/%s_%s.csv",outdir,baseNAME[2].c_str(),ext);
  fstream *fileout = new fstream(buffer,ios_base::out);

  (*fileout) << "N: " << N << endl;

  if( printALT ){
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "OR" << dels[0] << "95% CI <" << dels[0] << "95% CI > " << dels[0] << "p.value" << dels[0] << "p.adjusted" << dels[0] << "BC" << dels[0] << "p.value(ALT)" << dels[0] << "p.adjusted(ALT)" << dels[0] << "BC(ALT)" << endl;
  } else {
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "OR" << dels[0] << "95% CI <" << dels[0] << "95% CI > " << "p.value" << dels[0] << "p.adjusted" << dels[0] << "BC" << endl;
  }

  //linear indexing, size K is rows (A) x cols (B). 
  K=A*B;
  for(k=0; k<K; k++){

    a = floor(k/B);//row index
    b = k % B;     //col index

    string starsEO = "";
    string starsET = "";	
    string starsDO = "";
    string starsDT = "";
	
    for(i=0; i<SIGMASIZE; i++){
      if( p_values  [(a*B)+b] <= bonferroni[i] ){ starsEO += "*";  }
      if( p_valuesT [(a*B)+b] <= bonferroni[i] ){ starsET += "*";  }
      if( p_valuesD [(a*B)+b] <= bonferroni[i] ){ starsDO += "*";  }
      if( p_valuesDT[(a*B)+b] <= bonferroni[i] ){ starsDT += "*";  }
    }

    
    //--- Hypergeometric mean, or expected overlap
    double mn = 0.0;
    mn        = double(ANNOS[indexA][a].K) * double(ANNOS[indexB][b].K);
    mn       /= double(N);

    //--- Odds Ratio
    double ors = 0;
    double orL = 0;
    double orU = 0;

    double ORa = double(overlap[(a*B)+b]);
    double ORb = double(ANNOS[indexB][b].K - overlap[(a*B)+b]);
    double ORc = double(ANNOS[indexA][a].K - overlap[(a*B)+b]);
    double ORd = double(N - ANNOS[indexA][a].K + overlap[(a*B)+b] - ANNOS[indexB][b].K);
	
    calculateOddsRatio( ORa, ORb, ORc, ORd, ors, orL, orU );
    //---
    
    
    if( printALT ){
      (*fileout) <<  (printID ? ANNOS[indexA][a].annoID : ANNOS[indexA][a].annoDES) << dels[0] << ANNOS[indexA][a].K << dels[0] << (printID ? ANNOS[indexB][b].annoID : ANNOS[indexB][b].annoDES) << dels[0] << ANNOS[indexB][b].K << dels[0] << overlap[(a*B)+b] << dels[0] << mn << dels[0] << ors << dels[0] << orL << dels[0] << orU << dels[0] << (printTwoSided ? p_valuesT[(a*B)+b] : p_values[(a*B)+b]) << dels[0] << (printTwoSided ? padjustedT[(a*B)+b] : padjusted[(a*B)+b]) << dels[0] << (printTwoSided ? starsET : starsEO) << dels[0] << (printTwoSided ? p_valuesDT[(a*B)+b] : p_valuesD[(a*B)+b]) << dels[0] << (printTwoSided ? padjustedDT[(a*B)+b] : padjustedD[(a*B)+b]) << dels[0] << (printTwoSided ? starsDT : starsDO) << endl;
    } else {      
      (*fileout) <<  (printID ? ANNOS[indexA][a].annoID : ANNOS[indexA][a].annoDES) << dels[0] << ANNOS[indexA][a].K << dels[0] << (printID ? ANNOS[indexB][b].annoID : ANNOS[indexB][b].annoDES) << dels[0] << ANNOS[indexB][b].K << dels[0] << overlap[(a*B)+b] << dels[0] << mn << dels[0] << ors << dels[0] << orL << dels[0] << orU << dels[0] << (printTwoSided ? p_valuesT[(a*B)+b] : p_values[(a*B)+b]) << dels[0] << (printTwoSided ? padjustedT[(a*B)+b] : padjusted[(a*B)+b]) << dels[0] << (printTwoSided ? starsET : starsEO) << endl;
    }

  }
 
  (*fileout) << "" << endl;
  
  fileout->close();
      
}

void NetworkEnrichment::printOverlapinNetwork( int indexA, int indexB, int indexC, const char *outdir, const char *ext, bool printFDR){
  
  int i,k,K,a,b,c,A,B,C;

  A = Fsize[indexA];
  B = Fsize[indexB];
  C = Fsize[indexC]; 
  
  for(i=0; i<SIGMASIZE; i++){
    bonferroni[i] = sigma[i]/(double)(A*B*C);//Bonferroni at sigma[i]
  }

  
  char buffer  [BUFFERSIZE];
  sprintf(buffer ,"%s/%s_%s.csv",outdir,baseNAME[2].c_str(),ext);
  fstream *fileout = new fstream(buffer,ios_base::out);

  (*fileout) << "N: " << N << endl;

  if( printALT ){
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "annotation 3" << dels[0] << "n3" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "p.value" << dels[0] << "p.adjusted" << dels[0] << "BC" << dels[0] << "p.value(ALT)" << dels[0] << "p.adjusted(ALT)" << dels[0] << "BC(ALT)" << endl;
  } else {
    (*fileout) << "annotation 1" << dels[0] << "n1" << dels[0] << "annotation 2" << dels[0] << "n2" << dels[0] << "annotation 3" << dels[0] << "n3" << dels[0] << "actual overlap" << dels[0] << "expected overlap" << dels[0] << "p.value" << dels[0] << "p.adjusted" << dels[0] << "BC" << endl;
  }

  //linear indexing, size K is rows (A) x cols (B) x depth (C). 
  K=A*B*C;
  for(k=0; k<K; k++){

    c = k % C;                 //depth index
    b = ((k-c)/C) % B;         //col   index
    a = floor(((k-c)/C - b)/B);//row   index

    string starsEO = "";
    string starsET = "";	
    string starsDO = "";
    string starsDT = "";
	
    for(i=0; i<SIGMASIZE; i++){
      if( p_values  [c+C*(b+B*a)] <= bonferroni[i] ){ starsEO += "*";  }
      if( p_valuesT [c+C*(b+B*a)] <= bonferroni[i] ){ starsET += "*";  }
      if( p_valuesD [c+C*(b+B*a)] <= bonferroni[i] ){ starsDO += "*";  }
      if( p_valuesDT[c+C*(b+B*a)] <= bonferroni[i] ){ starsDT += "*";  }
    }
    //--- Hypergeometric mean, or expected overlap
    //--- See eqn (8) Alex T. Kalinka, The probablility of drawing intersections: extending the hypergeometric distribution,
    // arXiv:1305.0717v5, (2014).
    double mn = 0.0;
    mn        = double(ANNOS[indexA][a].K) * double(ANNOS[indexB][b].K) * double(ANNOS[indexC][c].K);
    mn       /= double(N*N);
    

    if( printALT ){
      (*fileout) <<  (printID ? ANNOS[indexA][a].annoID : ANNOS[indexA][a].annoDES) << dels[0] << ANNOS[indexA][a].K << dels[0] << (printID ? ANNOS[indexB][b].annoID : ANNOS[indexB][b].annoDES) << dels[0] << ANNOS[indexB][b].K << dels[0] << (printID ? ANNOS[indexC][c].annoID : ANNOS[indexC][c].annoDES) << dels[0] << ANNOS[indexC][c].K << dels[0] << overlap[c+C*(b+B*a)] << dels[0] << mn << dels[0] << (printTwoSided ? p_valuesT[c+C*(b+B*a)] : p_values[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? padjustedT[c+C*(b+B*a)] : padjusted[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? starsET : starsEO) << dels[0] << (printTwoSided ? p_valuesDT[c+C*(b+B*a)] : p_valuesD[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? padjustedDT[c+C*(b+B*a)] : padjustedD[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? starsDT : starsDO) << endl;  
    } else {      
      (*fileout) <<  (printID ? ANNOS[indexA][a].annoID : ANNOS[indexA][a].annoDES) << dels[0] << ANNOS[indexA][a].K << dels[0] << (printID ? ANNOS[indexB][b].annoID : ANNOS[indexB][b].annoDES) << dels[0] << ANNOS[indexB][b].K << dels[0] << (printID ? ANNOS[indexC][c].annoID : ANNOS[indexC][c].annoDES) << dels[0] << ANNOS[indexC][c].K << overlap[c+C*(b+B*a)] << dels[0] << dels[0] << mn << dels[0] << (printTwoSided ? p_valuesT [c+C*(b+B*a)] : p_values[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? padjustedT[c+C*(b+B*a)] : padjusted[c+C*(b+B*a)]) << dels[0] << (printTwoSided ? starsET : starsEO) << endl;
    }

  }
 
  (*fileout) << "" << endl;
  
  fileout->close();
      
}

//--------------------------------------------------------



//--------------------------------------------------------
//  HYPERGEOMETRIC TEST SETUP/DRIVER FUNCTIONS
//--------------------------------------------------------
int NetworkEnrichment::calculateOverlapinCommunities(bool runPermutations, const char* outdir, const char * ext, bool runFDR, bool alternativePrint, bool singlePerm ){  

  
  //N -> Number of genes, network size
  //M -> Number of communities
  //F -> Number of annotation types

  int i,j,k,p,m,f,K;

  //studies is a N (number of genes) x F (number annotation types) matrix
  cout << "Set memory for studies, " << N << "x" << F << "." << endl;
  studies = (double*)calloc(N*F,sizeof(double));  

  //linear indexing, size K is rows (N) x cols (F).
  K=N*F;
  //intialize numbers
  for(k=0; k<K; k++){ studies[k]=0; }

  //reset this first
  for(f=0; f<F; f++){ ANNOS[ANNOindex][f].K=0; }

  //find annotation types overlaping with genes in each cluster
  for( i=0; i<N; i++ ){
    int id = Clist[i].ID;
    for( k=0; k<Alines[ANNOindex]; k++ ){
      if( Alist[ANNOindex][k].ID == id ){
	for( j=0; j<F; j++ ){
	  if( strcmp(Alist[ANNOindex][k].annoID, ANNOS[ANNOindex][j].annoID) == 0 ){	  
	    studies[(i*F)+j]++; 
	    ANNOS[ANNOindex][j].K++; //tally of annotation type in network 
	  }
	}
      }
    }
  }

  //permute is a M (number of communities) x F (number annotation types) matrix
  cout << "Set memory for p_values & permute, " << M << "x" << F << "." << endl;
  p_values       = (double*)calloc(M*F,sizeof(double));
  padjusted      = (double*)calloc(M*F,sizeof(double));      
  permute        = (double*)calloc(M*F,sizeof(double));
  
  overlap        = (double*)calloc(M*F,sizeof(double)); 
  muab           = (double*)calloc(1,sizeof(double));//dummy

  p_valuesD      = (double*)calloc(M*F,sizeof(double));
  padjustedD     = (double*)calloc(M*F,sizeof(double));      
  permuteD       = (double*)calloc(M*F,sizeof(double));

  p_valuesDT     = (double*)calloc(M*F,sizeof(double));
  padjustedDT    = (double*)calloc(M*F,sizeof(double));      
  permuteDT      = (double*)calloc(M*F,sizeof(double));
  
  p_valuesT      = (double*)calloc(M*F,sizeof(double));
  padjustedT     = (double*)calloc(M*F,sizeof(double));      
  permuteT       = (double*)calloc(M*F,sizeof(double));
  
  //linear indexing, size K is rows (N) x cols (F).
  K=M*F;
  for(k=0; k<K; k++){
    p_values[k]   =0.0;
    padjusted[k]  =0.0;
    permute[k]    =0.0;

    overlap[k]    =0.0;

    p_valuesD[k]  =0.0;
    padjustedD[k] =0.0;
    permuteD[k]   =0.0;

    p_valuesDT[k] =0.0;
    padjustedDT[k]=0.0;
    permuteDT[k]  =0.0;
    
    p_valuesT[k]  =0.0;
    padjustedT[k] =0.0;
    permuteT[k]   =0.0;
  }
   
  comSIZE = (int*)calloc(M,sizeof(int));
  for(m=0; m<M; m++){
    comSIZE[m] = 0;
    comSIZE[m] = COMS[m].second;
  }

  annoSIZE = (int*)calloc(F,sizeof(int));
  for(f=0; f<F; f++){
    annoSIZE[f] = 0;
    annoSIZE[f] = ANNOS[ANNOindex][f].K;
  }

  geneCOM = (int*)calloc(N,sizeof(int));
  for(i=0; i<N; i++){
    geneCOM[i] = 0;
    geneCOM[i] = Clist[i].K;
  }

  //Overlap between genes in annotation set and communities
  overlapinNetwork();

  //Hypergeomertic test only
  cout << "calculate p-values...";
  overlapinComsHypergeometricTest();
  cout << "done." << endl;

  //FDR
  //if( runFDR ){
  //cout << "calculate adjust pvalues...";
  //calculateFDR();
  //cout << "done." << endl;
  //}

  if(runPermutations){//Hypergeomertic test & permutation study

    if( singlePerm ){

      cout << "run single permutation test...";
      
      permutation( 1.0 );//permute ids in the annotation list    
      overlapinComsHypergeometricTestRnd( false );//calculate & record permuated p-values

      cout << "done." << endl;      
      
    } else {
    
      cout << "permutation test, ints=" << NoP << "...";
      
      for(p=0; p<NoP; p++){      
	permutation( (double)(p+1) );//permute ids in the annotation list    
	overlapinComsHypergeometricTestRnd();//calculate & record permuated p-values
      }
  	
      cout << "done." << endl;

    }

  }      

  //FDR
  if( runFDR ){
    cout << "calculate adjust pvalues...";
    calculateFDR();
    cout << "done." << endl;
  }

  
  //print Hypergeomertic test & permutation study  
  printOverlapinCommunities( outdir, ext, runFDR, singlePerm );  
    
  if(alternativePrint){
    printOverlapinCommunitiesAlt( outdir, ext, runFDR );
  }
  

  freedMemory[0] = false;
  freeMemory();

  return 0;
  
}


//Overlap between two annotation files in a community/Bridging Region.
int NetworkEnrichment::calculateOverlapinCommunities( int indexA, int indexB, const char* outdir, const char* ext, bool runFDR ){

  //N -> Number of genes, network size
  //M -> Number of communities
  //A -> Number of annotation types A
  //B -> Number of annotation types B

  int m,k,a,b,A,B,K;
  
  //check we have at least two annotation sets and
  // indices A and B are valid
  if( (Alist.size() >= 2) && (indexA >= 0) && (indexA < Alist.size())
                          && (indexB >= 0) && (indexB < Alist.size())) {

    //indexA ==> ROW index
    //indexB ==> COL index
    A = Fsize[indexA];
    B = Fsize[indexB];

    //p_values is a A*B (number of annotation types A*B) x M (number of communities) matrix
    cout << "Set memory for p_value, " << A*B << "x" << M << "." << endl;
    p_values    = (double*)calloc(A*B*M,sizeof(double));
    padjusted   = (double*)calloc(A*B*M,sizeof(double));

    p_valuesD   = (double*)calloc(A*B*M,sizeof(double));
    padjustedD  = (double*)calloc(A*B*M,sizeof(double));

    p_valuesDT  = (double*)calloc(A*B*M,sizeof(double));
    padjustedDT = (double*)calloc(A*B*M,sizeof(double));

    p_valuesT   = (double*)calloc(A*B*M,sizeof(double));
    padjustedT  = (double*)calloc(A*B*M,sizeof(double));      

    //overlap of annotation types A and B in communities
    muCab       = (double*)calloc(A*B*M,sizeof(double));

    //product of annotation types A and B in communities
    nab         = (double*)calloc(A*B*M,sizeof(double));
    
    //overlap is a A (number of annotation types A) x B (number annotation types B) matrix
    cout << "Set memory for overlap, " << A << "x" << B << "." << endl;
    overlap = (double*)calloc(A*B,sizeof(double));
    muab    = (double*)calloc(A*B,sizeof(double));

    //set dummy value
    permute  = (double*)calloc(1,sizeof(double));
    annoSIZE = (int*)calloc(1,sizeof(int));
    geneCOM  = (int*)calloc(1,sizeof(int));    

    //linear indexing, size K is rows (A*B) x cols (M). 
    K=(A*B)*M;
    for(k=0; k<K; k++){
      p_values[k]    = 0;
      padjusted[k]   = 0;

      p_valuesD[k]   = 0;
      padjustedD[k]  = 0;

      p_valuesDT[k]  = 0;
      padjustedDT[k] = 0;

      p_valuesT[k]   = 0;
      padjustedT[k]  = 0;

      muCab[k]       = 0;
      nab[k]         = 0;
    }

    if( calRelDist ){
      //calculate the intersection distance
      p_dist      = (double*)calloc(A*B*M,sizeof(double));
      padjustedRD = (double*)calloc(A*B*M,sizeof(double));
      reldist     = (double*)calloc(A*B*M,sizeof(double));
      for(k=0; k<K; k++){
	p_dist[k]      = 0;
	padjustedRD[k] = 0;
	reldist[k]     = 0;
      }
    }
    
     //linear indexing, size K is rows (A) x cols (B). 
    K=A*B;
    for(k=0; k<K; k++){
      overlap[k] = 0;
      muab[k]    = 0;
    }   
    
    comSIZE = (int*)calloc(M,sizeof(int));
    for(m=0; m<M; m++){
      comSIZE[m] = 0;
      comSIZE[m] = COMS[m].second;
    }

    //Overlap between genes in two annotation sets
    overlapinNetork(indexA, indexB);
    
    cout << "calculate overlap p-values...";
    overlapinComsHypergeometricTest(indexA, indexB);
    cout << "done." << endl;

    //FDR
    if( runFDR ){
      cout << "calculate adjust pvalues...";
      calculateFDR(2, indexA, indexB);
      cout << "done." << endl;
    }

    //Print results
    printOverlapinCommunities(indexA, indexB, outdir, ext, runFDR );
    
  } else {
    cout << "indexA or indexB not valid." << endl;
    return 1;
  }

  freedMemory[0] = false;
  freeMemory();
  
  cout << "done." << endl;
  return 0;
  
}


//Overlap between two annotation files in the network
int NetworkEnrichment::calculateOverlapinNetwork( int indexA, int indexB, const char* outdir, const char* ext, bool runFDR ){


  //N -> Number of genes, network size
  //M -> Number of communities
  //A -> Number of annotation types A
  //B -> Number of annotation types B

  int k,K,A,B;
  
  //check we have at least two annotation sets and
  // indices A and B are valid
  if( (Alist.size() >= 2) && (indexA >= 0) && (indexA < Alist.size())
                          && (indexB >= 0) && (indexB < Alist.size())) {

    //indexA ==> ROW index
    //indexB ==> COL index
    A = Fsize[indexA];
    B = Fsize[indexB];

    //p_values is a A (number of annotation types A) x B (number annotation types B) matrix 
    cout << "Set memory for p_value, " << A << "x" << B << "." << endl;
    p_values    = (double*)calloc(A*B,sizeof(double));
    padjusted   = (double*)calloc(A*B,sizeof(double));

    p_valuesD    = (double*)calloc(A*B,sizeof(double));
    padjustedD   = (double*)calloc(A*B,sizeof(double));

    p_valuesDT   = (double*)calloc(A*B,sizeof(double));
    padjustedDT  = (double*)calloc(A*B,sizeof(double));

    p_valuesT    = (double*)calloc(A*B,sizeof(double));
    padjustedT   = (double*)calloc(A*B,sizeof(double));      
    
    //overlap is a A (number of annotation types A) x B (number annotation types B) matrix
    cout << "Set memory for overlap, " << A << "x" << B << "." << endl;
    overlap = (double*)calloc(A*B,sizeof(double));
    muab    = (double*)calloc(A*B,sizeof(double));  
    
    //set dummy value
    permute  = (double*)calloc(1,sizeof(double));
    comSIZE  = (int*)calloc(1,sizeof(int));
    annoSIZE = (int*)calloc(1,sizeof(int));
    geneCOM  = (int*)calloc(1,sizeof(int));

    //linear indexing, size K is rows (A) x cols (B). 
    K=A*B;
    for(k=0; k<K; k++){
      p_values[k]    = 0;
      padjusted[k]   = 0;

      p_valuesD[k]   = 0;
      padjustedD[k]  = 0;

      p_valuesDT[k]  = 0;
      padjustedDT[k] = 0;

      p_valuesT[k]   = 0;
      padjustedT[k]  = 0;

      overlap[k]     = 0;
      muab[k]        = 0;
    }

    overlapinNetork(indexA, indexB);
    
    cout << "calculate overlap p-values...";
    overlapinNetHypergeometricTest(indexA, indexB);
    cout << "done." << endl;

    //FDR
    if( runFDR ){
      cout << "calculate adjust pvalues...";
      calculateFDR(3, indexA, indexB);
      cout << "done." << endl;
    }
    
    printOverlapinNetwork(indexA, indexB, outdir, ext, runFDR);
    
  } else {
    cout << "indexA or indexB not valid." << endl;
    return 1;
  }

  freedMemory[0] = false;
  freeMemory();
  
  cout << "done." << endl;
  return 0;
   

}

//Overlap between three annotation files in the network
int NetworkEnrichment::calculateOverlapinNetwork( int indexA, int indexB, int indexC, const char* outdir, const char* ext, bool runFDR ){


  //N -> Number of genes, network size
  //M -> Number of communities
  //A -> Number of annotation types A
  //B -> Number of annotation types B
  //C -> Number of annotation types C

  int k,K,a,b,c,A,B,C;

  
  //check we have at least three annotation sets and
  // indices A, B and C are valid
  if( (Alist.size() >= 3) && (indexA >= 0) && (indexA < Alist.size())
                          && (indexB >= 0) && (indexB < Alist.size()) 
                          && (indexC >= 0) && (indexC < Alist.size())){
    
    //indexA ==> ROW index
    //indexB ==> COL index
    //indexC ==> COL index
    A = Fsize[indexA];
    B = Fsize[indexB];
    C = Fsize[indexC];

    //p_values is a A (number of annotation types A) x B (number annotation types B) x C matrix 
    cout << "Set memory for p_value, " << A << "x" << B << "x" << C << "." << endl;
    p_values     = (double*)calloc(A*B*C,sizeof(double));
    padjusted    = (double*)calloc(A*B*C,sizeof(double));

    p_valuesD    = (double*)calloc(A*B*C,sizeof(double));
    padjustedD   = (double*)calloc(A*B*C,sizeof(double));

    p_valuesDT   = (double*)calloc(A*B*C,sizeof(double));
    padjustedDT  = (double*)calloc(A*B*C,sizeof(double));

    p_valuesT    = (double*)calloc(A*B*C,sizeof(double));
    padjustedT   = (double*)calloc(A*B*C,sizeof(double));
    
    //overlap is a A (number of annotation types A) x B (number annotation types B) X C matrix
    cout << "Set memory for overlap, " << A << "x" << B << "x" << C << "." << endl;
    overlap = (double*)calloc(A*B*C,sizeof(double));
    muab    = (double*)calloc(A*B*C,sizeof(double));  
    
    //set dummy value
    permute  = (double*)calloc(1,sizeof(double));
    comSIZE  = (int*)calloc(1,sizeof(int));
    annoSIZE = (int*)calloc(1,sizeof(int));
    geneCOM  = (int*)calloc(1,sizeof(int));

    //linear indexing, size K is rows (A) x cols (B) x depth (C). 
    K=A*B*C;
    for(k=0; k<K; k++){
      p_values[k]    = 0;
      padjusted[k]   = 0;

      p_valuesD[k]   = 0;
      padjustedD[k]  = 0;

      p_valuesDT[k]  = 0;
      padjustedDT[k] = 0;

      p_valuesT[k]   = 0;
      padjustedT[k]  = 0;

      overlap[k]     = 0;
      muab[k]        = 0;
    }
    
    overlapinNetork(indexA, indexB, indexC);
    
    cout << "calculate overlap p-values...";
    overlapinNetHypergeometricTest(indexA, indexB, indexC);
    cout << "done." << endl;

    //FDR    
    if( runFDR ){
      cout << "calculate adjust pvalues...";
      calculateFDR(4, indexA, indexB, indexC);
      cout << "done." << endl;
    }
    

    printOverlapinNetwork(indexA, indexB, indexC, outdir, ext, runFDR);
    
  } else {
    cout << "indexA or indexB not valid." << endl;
    return 1;
  }

  freedMemory[0] = false;
  freeMemory();
  
  
  cout << "done." << endl;
  return 0;
   

}


//Overlap between annotation file and communities in the network
int NetworkEnrichment::calculateOverlapinCommunities( int Index, const char* outdir, const char* ext, bool runFDR ){

  int cal;
  
  setANNOindex( Index );

  cal = calculateOverlapinCommunities(false, outdir, ext, runFDR, true);

  return cal;
  
}



