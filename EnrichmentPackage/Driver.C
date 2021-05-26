//----------------------------------------------
// Package Name       : EnrichmentPackage  
// Package Version    : v3r1
// Package Description: Use of the Hypergeometric distribution  
//                    : to calculate enrichment in clustered PPI networks
// Date               : 2016
// Author             : Colin D Mclean <Colin.D.Mclean@ed.ac.uk>
// Copyright (C) 2016 Colin Mclean 
//----------------------------------------------
//      Package Description
//----------------------------------------------
// Use of the Hypergeometric distribution to calculate enrichment
// in clustered PPI networks, built upon:
// [1] the gamma function (and class) as given in Numerical Recipies 3rd Edition W. Press, S. Teukolsky, W. Vetterling, B. Flannery.
// [2] M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078.
// [3] Pocklington A, Cumiskey D, Armstrong D, Grant S: The proteomes of neurotransmitter receptor complexes from modular networks
//     with distributed functionality underlying plasticity and behaviour, MSB, 2, (2006).
// [4] Alex T. Kalinka, The probablility of drawing intersections: extending the hypergeometric distribution, arXiv:1305.0717v5, (2014).
// [5] Benjamini, Y., and Hochberg, Y. Controlling the false discovery rate:  a practical and powerful approach to multiple testing.
//     Journal of the Royal Statistical Society Series B 57 (1995), 289–300.
// [6] Benjamini, Y., and Liu, W. A step-down multiple hypotheses testing procedure that controls the false discovery rate under independence.
//     Journal of Statistical Planning and Inference 82 (1999), 163–170.
// [7] Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188.
//----------------------------------------------
//      GNU General Public Licenses v3 
//----------------------------------------------
// This program is free software: you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the 
// Free Software Foundation, either version 3 of the License, or (at your 
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// (GNU_GPL_v3)  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------
//      Funding Acknowledgement
//----------------------------------------------
// This open source software code was developed in part or in whole in the Human Brain Project, funded from the European Union’s Horizon 2020
// Framework Programme for Research and Innovation under the Specific Grant Agreement No. 720270 (Human Brain Project SGA1).
//-----------------------------------------------
/////////////////////////////////////////////////
//-----------------------------------------------
#include "Headers.h"
#include "NetworkEnrichment.h"

using namespace std;


//--- METHOD DECLARATIONS
//-------------------------------------------------------------------------------------
int readDirectory    ( const char *, vector<string> &);
void printLicence    ( const char * );
void printHelpMessage( const char * );
bool matchUsersInput ( const char *, vector<string> );
void setupMatrices ( );
void freeMemory();

//--- GLOBAL VARIABLES
//-------------------------------------------------------------------------------------



enum {NONE};



//--- MAIN PROGRAM
//-------------------------------------------------------------------------------------
int main(int argc, char * argv[]) {

  int i, seed, nop, oo1, cal1, cal2, cal3, cal4, Option;
  double pesudoCount;
  string if_help;
  const char *PERM=NULL;
  const char *ADDoffset=NULL;
  const char *FDRmeth=NULL;
  const char *PRINTid=NULL;
  const char *Ext=NULL;
  const char *file_names;
  const char *OUTDIR=NULL;
  vector<string> SetName;
  vector<string> INfiles;
  vector<int>    MINov;

  seed            = -1;
  nop             = -1;

  pesudoCount     = -1.0;
  
  bool useFDR       = false;
  bool usePerm      = true;
  bool usePrintID   = false;
  bool usePrintAn   = false;
  bool usePrintCnew = false;
  
  //bool usePrintFC     = false;
  //bool noPrintMeanMu  = false;
  bool noPrintALT     = false;
  
  bool useOneSided = false;
  bool useTwoSided = false;
  //bool useRelDist  = false;
  bool useMaxSS    = false;

  bool useRCfisher = false;
  bool useChi2     = false;  
  
  //Set possible 'positive' user choices
  vector<string> POSchoices;
  POSchoices.push_back("");
  POSchoices.push_back("y");
  POSchoices.push_back("yes");
  POSchoices.push_back("t");
  POSchoices.push_back("true");

  //Set possible 'negative' user choices
  vector<string> NEGchoices;
  NEGchoices.push_back("");
  NEGchoices.push_back("n");
  NEGchoices.push_back("no");
  NEGchoices.push_back("f");
  NEGchoices.push_back("false"); 

  //Set possible FDR methods
  vector<string> FDRchoices;
  FDRchoices.push_back("");
  FDRchoices.push_back("BH");
  FDRchoices.push_back("BY");
  FDRchoices.push_back("BL");

  
  
  if ( argc == 1 ){ 
    printHelpMessage( argv[0] );
  }

  Option = 1;
  int a  = 1;
  while( a <= argc -1 ){

    if( strcmp(argv[a],"-h") == 0 || strcmp(argv[a],"-help") == 0 ){
     printHelpMessage( argv[0] );
    } else if( strcmp(argv[a],"-noPerm") == 0 ){
      //PERM = argv[++a];
      usePerm = false;
      //} else if( strcmp(argv[a],"-offset") == 0 ){
      //ADDoffset = argv[++a];
    } else if( strcmp(argv[a],"-setFDR") == 0 ){
      FDRmeth = argv[++a];
    } else if( strcmp(argv[a],"-setSEED") == 0 ){
      seed = (int)atoi(argv[++a]);
    } else if( strcmp(argv[a],"-setNoP") == 0 ){
      nop = (int)atoi(argv[++a]);
    } else if( strcmp(argv[a],"-setPesudoC") == 0 ){
      pesudoCount = (double)atof(argv[++a]);
    } else if( strcmp(argv[a],"-printID") == 0 ){
      //PRINTid = argv[++a];
      usePrintID = true;
    } else if( strcmp(argv[a],"-printAn") == 0 ){
      usePrintAn = true;
    } else if( strcmp(argv[a],"-printCnew") == 0 ){
      usePrintCnew = true;
      //} else if( strcmp(argv[a],"-printFC") == 0 ){
      //usePrintFC = true;
      //} else if( strcmp(argv[a],"-noExMu") == 0 ){
      //noPrintMeanMu = true;
    } else if( strcmp(argv[a],"-noALT") == 0 ){
      noPrintALT = true;
    } else if( strcmp(argv[a],"-onesided") == 0 ){
      useOneSided = true;
    } else if( strcmp(argv[a],"-twosided") == 0 ){
      useTwoSided = true;
    } else if( strcmp(argv[a],"-maxSS") == 0 ){
      useMaxSS  = true;
      //} else if( strcmp(argv[a],"-reldist") == 0 ){
      //useRelDist = true;
    } else if( strcmp(argv[a],"-RCfisher") == 0 ){
      useRCfisher = true;
    } else if( strcmp(argv[a],"-Chi2") == 0 ){
      useChi2 = true;
    } else if( strcmp(argv[a],"-setN") == 0 ){
      SetName.push_back(argv[++a]);
    } else if( strcmp(argv[a],"-opt") == 0 ){
      Option = (int)atoi(argv[++a]);
    } else if( strcmp(argv[a],"-minOV") == 0 ){
      MINov.push_back((int)atoi(argv[++a]));
    } else if( strcmp(argv[a],"-ext") == 0 ){
      Ext = argv[++a];
    } else if( strcmp(argv[a],"-Comfile") == 0 ){
      INfiles.push_back(argv[++a]);
    } else if( strcmp(argv[a],"-Annofile") == 0 ){
      INfiles.push_back(argv[++a]);    
    } else if( strcmp(argv[a],"-o") == 0 ){
      OUTDIR = argv[++a];
    } else if( strcmp(argv[a],"-lic") == 0 ){
      printLicence( argv[0] );
    }
    ++a;
  }

  
  if( Ext == NULL ){ Ext = "TEMP"; }
  
  //create the output directory
  if( OUTDIR == NULL ){ OUTDIR = "OUT"; }
  oo1 = fileReader::writeDirectory( OUTDIR );

  if( oo1 != 0 ){ cout << "Directory " << OUTDIR << " exits. " << endl; }
  else { cout << "Creating directory " << OUTDIR << endl; } 

  //Create Network Enrichment object and pass our input files to it.
  NetworkEnrichment *enrD = new NetworkEnrichment( INfiles );

  //Before any calculations, check if we should off-set the community numbers.
  //For example, you'd invoke this if you're community numbers started from '0'.
  if( matchUsersInput(ADDoffset, POSchoices) ){ enrD->setKOffset(); }  

  //
  //if( noPrintMeanMu ){ enrD->setFoldChange( false ); }

  //
  if( noPrintALT ){ enrD->setALT( false ); }

  //
  //if( useRelDist ){ enrD->setRelDist( true ); }

  //
  if( useRCfisher ){ enrD->setRCfisher( true ); }  

  //
  if( useChi2 ){ enrD->setChi2( true ); }  
  
  //
  //if( usePrintFC ){ enrD->setFoldChange( true ); }    

  //Whether to print annotation IDs, or description
  if( usePrintCnew ){ enrD->setPrintCNEW( true ); }    
  
  //Whether to print annotation IDs, or description
  if( usePrintID ){ enrD->setPrintID( true ); }    

  //Whether to print annotation Size
  if( usePrintAn ){ enrD->setPrintAn( true ); }  

  //calculate depletion
  if( useOneSided ){ enrD->oneSided(); }

  //calculate two-sided
  if( useTwoSided ){ enrD->twoSided(); }

   //calculate two-sided
  if( useMaxSS ){ enrD->maxSS(); }
  
  //Check if the user what's permutations
  //if( matchUsersInput(PERM, NEGchoices) ){ usePerm = false; }  
  
  //Set which FDR method to use: BY (default), BH, BL. 
  if( matchUsersInput(FDRmeth,FDRchoices) ){ enrD->setFDRmethod( FDRmeth ); useFDR = true; }

  for(i=0; i<MINov.size(); i++){ enrD->setMINOVERLAP (i, MINov[i]); }
  MINov.clear();  
  
  //
  if( Option == 1 ){
    if( seed != -1        ){ enrD->seedOffSet(true, seed);        }
    if( nop  != -1        ){ enrD->setNoP( nop );                 }
    if( pesudoCount != -1 ){ enrD->setPesudoCount( pesudoCount ); }
    
    cal1 = enrD->calculateOverlapinCommunities(usePerm, OUTDIR, Ext, useFDR, false, false);
  }

  if( Option == 2 ){
    cal2 = enrD->calculateOverlapinCommunities(1, 0, OUTDIR, Ext);
  }
  
  if( Option == 3 ){

    if( SetName.size() == 0 ){ SetName.push_back("TEMPA"); SetName.push_back("TEMPB"); }
    
    cal3 = enrD->calculateOverlapinNetwork(1, 0,   OUTDIR, Ext);
    //cal3 = enrD->calculateOverlapinCommunities(0, OUTDIR, SetName[0].c_str());
    //cal3 = enrD->calculateOverlapinCommunities(1, OUTDIR, SetName[1].c_str());
  }
  

   if( Option == 4 ){

     if( SetName.size() == 0 ){ SetName.push_back("TEMPA"); SetName.push_back("TEMPB"); SetName.push_back("TEMPC"); }
     
     cal4 = enrD->calculateOverlapinNetwork(0, 1, 2, OUTDIR, Ext);
     cal4 = enrD->calculateOverlapinCommunities(0,  OUTDIR, SetName[0].c_str());
     cal4 = enrD->calculateOverlapinCommunities(1,  OUTDIR, SetName[1].c_str());
     cal4 = enrD->calculateOverlapinCommunities(2,  OUTDIR, SetName[2].c_str());

   }

   //
   if( Option == 5 ){
     if( seed != -1 ){ enrD->seedOffSet(true, seed); }
     cal1 = enrD->calculateOverlapinCommunities(true, OUTDIR, Ext, useFDR, false, true);
   }

  
  delete enrD;  

  freeMemory();
  
  exit(1);

}

void printLicence( const char *arg1){

  cout << "***-----***" << endl;
  cout << "Network Enrichment Package (v3r1) " << " Copyright (C) " << " 2016 " << " Colin D Mclean " << endl;
  cout << "***-----***" << endl;
  cout << "This program comes with ABSOLUTELY NO WARRANTY;" << endl;
  cout << "This is free software, and you are welcome to redistribute it  under certain conditions; see file \"GNU_GPL_v3\" for details. " << endl;
  cout << "-----" << endl;
  cout << "THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT  PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE  PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION. " << endl;
   cout << "***-----***" << endl;
   cout << " Funding Acknowledgement " << endl;
   cout << "***-----***" << endl;
   cout << " This open source software code was developed in part or in whole in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under the Specific Grant Agreement No. 720270 (Human Brain Project SGA1)." << endl;   
   cout << "-----" << endl;
   
  exit(1);

}

void printHelpMessage( const char *arg1){

  cout << arg1 << " arguments: " << endl;  
  cout << "-Comfile      : Community File for enrichment studies " << endl;
  cout << "-Annofile     : Annotation File for enrichment studies " << endl;
  cout << "-opt          : " << endl;
  cout << "              : [1] Run cluster enrichment on an annotation set" << endl;
  cout << "              : 2 = Run cluster enrichment on two annotation set sets" << endl;
  cout << "              : 3 = Run network enrichment on two annotation set sets" << endl;
  cout << "              : 4 = Run network enrichment on three annotation set sets" << endl;
  cout << "-setFDR       : Set FDR method to adjust p-values" << endl;
  cout << "              : [BY] Default is to run Benjamini and Yekutieli FDR method" << endl;
  cout << "              : BH = run Benjamini and Hochberg FDR method" << endl;
  cout << "              : BL = run Benjamini and Liu FDR method" << endl;
  cout << "-setSEED      : add value to our random number seed" << endl;
  cout << "              : [1 (int)] Default" << endl;
  cout << "-setPesudoC   : avoid zero p-values in the permutation test" << endl;
  cout << "              : [1.0 (double)] Default" << endl;
  cout << "-setNoP       : number of permutation tests" << endl;
  cout << "              : [1000 (int)] Default" << endl;  
  cout << "-minOV        : Set minimum values for: " << endl;
  cout << "              : 1st call = set min community size " << endl;
  cout << "              : 2nd call = set min annotation size " << endl;
  cout << "              : 3rd call = set min annotation overlap size in a community " << endl;
  //cout << "-offset       : [y,yes] invoke if community in Comfile start from 0" << endl;
  cout << "-printCnew    : print new community numbers. " << endl;
  cout << "              : Internally we will rescale community numbers from 1 to (max(C) - min(C))." << endl;
  cout << "              : Default is to print original community numbers in -Comfile." << endl;
  cout << "-printID      : to print annotation IDs. " << endl;
  cout << "              : Default is to print annotation description." << endl;
  cout << "-printAn      : to print annotation Size. " << endl;  
  cout << "              : Default is not too." << endl;
  //cout << "-printFC      : print fold-change. " << endl;  
  //cout << "              : Default is not too." << endl;
  //cout << "-noExMu       : don't print expected overlap. " << endl;  
  //cout << "              : Default is too." << endl;
  cout << "-noALT        : don't print alternative enrichment/depletion values. " << endl;  
  cout << "              : Default is too." << endl;
  cout << "-onesided     : print erichment / depletion using one-sided Fisher test. "  << endl;
  cout << "-twosided     : print enrichment / depletion using two-sided Fisher test. " << endl;
  cout << "              : Default is to calculate enrichment / depletion using two-sided Fisher test." << endl;
  cout << "-maxSS        :  When calculating the overlap between two anntation sets in communities, do we compare this against the 'maximum' sample space or 'standard' sample space." << endl;
  cout << "              : Default is to compare against 'standard' sample space." << endl;
  cout << "-RCfisher     : Calculating the exact Fisher test p.value on overlap between two anntation sets in communities, using R's RxC contingency table code." << endl;
  cout << "              : Default is not to use this." << endl;
  cout << "-Chi2         : Calculating the Chi2 test p.value on overlap between two anntation sets in communities, using the RxC contingency tables." << endl;
  cout << "              : Default is not to use this." << endl;
  /*cout << "-reldist      : calculate & print the relative distance between two annotation hypergeometric distributions. " << endl;
    cout << "              : Default is not to calculate & print the relative distance." << endl;*/
  cout << "-noPerm       : Don't calculate permutation on p.values. " << endl;
  cout << "              : Default is to calculate p.values permutations." << endl;
  cout << "-ext          : Add a description to end of output files " << endl;
  cout << "-o            : Output directory " << endl;
  cout << "-lic          : Print program license details and funding acknowledgement" << endl;  
  cout << "-----" << endl;
  cout << "EXAMPLE 1 : Run Cluster enrichment over one annotation set: e.g. disease set" << endl;
  cout << "-----" << endl;
  cout << "./run -opt 1 -Comfile ../Clustering/PPI_Presynaptic_Published/Spectral_communities.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv" << endl;
  cout << "-----" << endl;
  cout << "EXAMPLE 2 : Run Cluster enrichment over two annotation sets: e.g. disease and synaptic functional annotation sets" << endl;
  cout << "-----" << endl;
  cout << "./run -opt 2 -Comfile  ../Clustering/PPI_Presynaptic_Published/Spectral_communities.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv -Annofile ../Annotations/SCH_flatfile.csv" << endl;
  cout << "-----" << endl;
  cout << "EXAMPLE 3 : Run Network enrichment over three annotation sets: e.g. disease, synaptic function and cell type annotation sets" << endl;
  cout << "-----" << endl;
  cout << "./run -opt 4 -Comfile  ../Clustering/PPI_Presynaptic_Published/Spectral_communities.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv -Annofile ../Annotations/SCH_flatfile.csv -Annofile ../Annotations/celltypes_PMID27991900_L2.csv" << endl;
  cout << "-----" << endl;
  cout << "EXAMPLE 4 : Run a single Permutation test for cluster enrichment on annotation set" << endl;
  cout << "-----" << endl;
  cout << "./run -opt 5 -Comfile ../Clustering/PPI_Presynaptic_Published/Spectral_communities.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv" << endl;
  cout << "-----" << endl;
    cout << "-----" << endl;
  cout << "EXAMPLE 5 :  Run Network enrichment over two annotation sets: e.g. disease, synaptic function annotation sets" << endl;
  cout << "-----" << endl;
  cout << "./run -opt 3 -Comfile ../Clustering/PPI_Presynaptic_Published/Spectral_communities.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv -Annofile ../Annotations/SCH_flatfile.csv" << endl;
  cout << "-----" << endl;
 
  exit(1);

}

bool matchUsersInput( const char * TEST, vector<string> CHOICES ){

  int i;
  
  if( TEST == NULL ){ return false; }

  for(i=0; i<CHOICES.size(); i++){
    if( strcasecmp(TEST,CHOICES[i].c_str()) == 0 ){ return true; }
  }

  return false;
  
}
     
void setupMatrices(){}


void freeMemory(){}

