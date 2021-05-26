#include "Headers.h"
#include "fileReader.h"

//define guards, so headers are declare only once.
#ifndef BUILDSETS_H
#define BUILDSETS_H

//define useful list struct objects needed for annotation files
struct LISTst {
    char *annoID;
    char *annoDES;
    int   ID;
    int   K;
  };

struct sortpairDoubInt{
  bool operator()(const std::pair<double,int> &l, const std::pair<double,int> &r) const {
    return l.first < r.first;
    }
  };


struct sortpairIntInt{
  bool operator()(const std::pair<int,int> &l, const std::pair<int,int> &r) const {

    if( l.first < r.first ) return true;
    if( l.first > r.first ) return false;
    return l.second < r.second;    
    //return l.first < r.first;

    }
  };


struct sortpairStrInt{
  bool operator()(const std::pair<string,int> &l, const std::pair<string,int> &r) const {

    if( l.first.compare(r.first) < 0 ) return true; 
    if( l.first.compare(r.first) > 0 ) return false;
    return l.second < r.second;
    //return l.first.compare(r.first)<0;

    }
  };


 //Print _Nr x _Nc matrix _M
  static inline void printM( double *_M, int _Nr, int _Nc, const char* _Name ){

    int i,j,k,KK;
    
    //char buffer[250];
    //sprintf(buffer ,"Printing Matrix %s:",_Name);
    cout << "Printing Matrix: " << _Name << endl;

    KK=_Nr * _Nc;
    for(k=0; k<KK; k++){
    i = floor(k/_Nr);
    j = k % _Nr;
    if( _M[(i*_Nr)+j] != 0 ) cout << " " << _M[(i*_Nr)+j] << " ";
    else                     cout << " . ";

    if( j == (_Nr-1) )       cout << "" << endl;
    }
 
 };

class buildSets {

 public:
   buildSets();
   //void addSets(const char *[], int );
   void addSets(vector<string>, int );
  ~buildSets();

  void readAnnotationFile (const char *); 
 
 private:
  void    assignSpace();
  void    freeSpace();
  LISTst *freqofAnnolist( LISTst*, int, int & );
  void    freqofComslist( bool = false, int = 0 );
  LISTst *removeDuplicateIDs( LISTst *, int &, LISTst *, int );
  
  LISTst *createList( fileReader *, int, int );


  //data files
  int        Nfiles;  
  bool       freedfiles;
  vector<fileReader*> files;

  //annotation sets
  bool   freedClist;  //temp variables for community file  
  bool   freedAlist;  //temp variables for annotation files  
  bool   freedANNOS;

  //set file delineator(s)
  static const int DELSIZE = 1;
  char dels[DELSIZE];

  //set file header delineator(s)
  static const int HEADDELSIZE = 1;
  char HEADdel[HEADDELSIZE];
  
 protected:

  void printList( LISTst *, int );
  int  getKMin();
  void addKOffset( int = 0 );
  
  typedef std::pair<double,int> pairDoubInt;
  typedef std::pair<int,int>    pairIntInt;
  typedef std::pair<string,int> pairStrInt;

  typedef std::tuple<int,int,int> tripleInt;
  
  //for community file
  int     Mmin_old;
  int     Mmax_old;
  int     M_old;
  int     Mmin;
  int     Mmax;
  int     Clines;
  int     Ccols;
  LISTst *Clist;
  vector<tripleInt> COMS; 

  //for annotation files
  vector<int>     Alines;
  int             Acols;
  vector<LISTst*> Alist;
  vector<LISTst*> ANNOS;

  //make these protected
  int N;
  int M;
  int F;
  vector<int> Fsize;

};

#endif
