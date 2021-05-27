#include "buildSets.h"

buildSets::buildSets() {

  int i;
  
  this->Nfiles = 0;
  this->Clist  = 0;
  this->Clines = 0;
  this->Ccols  = 0;
  this->Acols  = 0;
   
  this->N        = 0;
  this->M        = 0;
  this->M_old    = 0;
  this->Mmin_old = 0;
  this->Mmax_old = 0;
  this->Mmin     = 0;
  this->Mmax     = 0;  
  this->F        = 0;
 
  
  this->freedfiles = false;
  this->freedClist = false;
  this->freedAlist = false;
  this->freedANNOS = false;

  this->dels[0]='\t';
  //this->dels[1]=',';
  //this->dels[2]=' ';

  this->HEADdel[0]='#';
  
}

void buildSets::addSets( const char *INfiles[], int Nfiles ) {
 
  int i,j,k,nlines,Ff;
  
  this->N        = 0;
  this->M        = 0;
  this->M_old    = 0;
  this->Mmin_old = 0;
  this->Mmax_old = 0;
  this->Mmin     = 0;
  this->Mmax     = 0;
  this->F        = 0;
 
  this->Nfiles     = Nfiles;
  this->freedfiles = false;
  this->freedClist = false;
  this->freedAlist = false;
  this->freedANNOS = false;

  this->dels[0]='\t';
  //this->dels[1]=',';
  //this->dels[2]=' ';

  this->HEADdel[0]='#';
  
  //create space for input files
  //cout << "Nfiles: " << Nfiles << endl;
  for( i=0; i<Nfiles; i++ ){
    cout << "reading " << INfiles[i] << endl;
    this->files.push_back(new fileReader((FILE*)fopen(INfiles[i],"r")));
  }
 
  //fill space with input files
  //fill community file
  //files[0] will/should always be the community file
  files[0]->fill_buffer();
  Clines = files[0]->get_nlines();
  Ccols  = 2;  
  Clist = createList( files[0], Clines, Ccols );

  //files[1,Nfiles] the annotation files
  Acols  = 3; 
  for( i=1; i<Nfiles; i++){  
    files[i]->fill_buffer();
    nlines = files[i]->get_nlines();
    this->Alines.push_back(nlines);
    this->Alist.push_back(createList( files[i], nlines, Acols ));
  }
  
  //can afford to delete input files here 
  if( freedfiles == false ){
    for( i=0; i<files.size(); i++){ 
      delete files[i]; vector<fileReader*>().swap(files); freedfiles=true; }
  }    

  //calculate unqiue number of annotations types & the sizes.
  freqofComslist();
  for( i=0; i<(Nfiles-1); i++ ){
    Ff = 0;
    this->ANNOS.push_back( freqofAnnolist(Alist[i], Alines[i], Ff) );
    this->Fsize.push_back( Ff );
  }

  //check for duplicate Entrez IDs in each annotation type, which may exist
  //in the input annotation files.
  for( i=0; i<(Nfiles-1); i++ ){
    Alist[i] = removeDuplicateIDs( Alist[i], Alines[i], ANNOS[i], Fsize[i] );    
  }

  N = Clines;
  M = COMS.size();
  F = Fsize[0];
    
}


void buildSets::readAnnotationFile(const char *annoFile ){
  
  int i, nlines, Ff;

  //reset annotation space
  freeSpace();
  freedAlist = false;
  freedfiles = false;

  //create space for annotation file
  cout << "reading " << annoFile << endl;
  this->files.push_back(new fileReader((FILE*)fopen(annoFile,"r")));
  
  files[0]->fill_buffer();
  nlines = files[0]->get_nlines();
  this->Alines.push_back(nlines);
  Acols  = 3; 
  this->Alist.push_back( createList(files[0], nlines, Acols) );

  //can afford to delete input files here 
  if( freedfiles == false ){
    for( i=0; i<files.size(); i++){ 
      delete files[i]; vector<fileReader*>().swap(files); freedfiles=true; }
  }    

  //calculate unqiue number of annotations types & the sizes.
  Ff = 0;
  this->ANNOS.push_back( freqofAnnolist(Alist[0], Alines[0], Ff) );
  this->Fsize.push_back(Ff);

  //check for duplicate Entrez IDs in each annotation type, which may exist
  //in the input annotation files.
  removeDuplicateIDs( Alist[0], Alines[0], ANNOS[0], Fsize[0] );    
    
}



LISTst* buildSets::createList( fileReader *infile, int NROWS, int NCOLS ){

  int i,j,k,start,end,length,row,col;
  vector<int>    del_pos;
  vector<string> tokens;

  char line[1000];
  char token[1000];
  char *nonspace;

  LISTst *LIST = (LISTst*)calloc(NROWS,sizeof(LISTst));

  //reset buffer 
  infile->reset_buffer();

  //read each line in file buffer
  row=0; col=0;
  while( infile->next_line(line) == 0 ){

    //skip lines starting with '#'
    if(line[0]==HEADdel[0]) continue;
  
    start  = 0;
    end    = 0;
    length = strlen(line);
    
    fileReader::findDelimiterPos(line, length, dels, DELSIZE, del_pos );

    //if no delemiters skip
    if(del_pos.size() != 0){

      for(k=0; k<del_pos.size(); k++){

	end=del_pos[k];
	
	fileReader::subString(line,token,(start+1),end);	  
	for(nonspace=token; *nonspace==' '; nonspace++);
	tokens.push_back(nonspace);
	start=end+1;//+1 to miss out the delimiter
      }
      
      if( start < length ){
	fileReader::subString(line,token,(start+1),length);
	for(nonspace=token; *nonspace==' '; nonspace++);
	tokens.push_back(nonspace);
      }
      
    }

    
    if(tokens.size() == 2){//community file
      LIST[row].ID = atoi(tokens[0].c_str());
      LIST[row].K = atoi(tokens[1].c_str());
    }

    if(tokens.size() == 3){//annotation file
      
      LIST[row].annoID = (char*)malloc((strlen(tokens[0].c_str())+1)*sizeof(char));
      strcpy(LIST[row].annoID,tokens[0].c_str());

      LIST[row].annoDES = (char*)malloc((strlen(tokens[1].c_str())+1)*sizeof(char));
      strcpy(LIST[row].annoDES,tokens[1].c_str());

      LIST[row].ID = atoi(tokens[2].c_str());
    }
    

    tokens.clear();
    del_pos.clear();
    col=0;
    row++;

  }

  vector<int>().swap(del_pos);  //free space of vector
  vector<string>().swap(tokens);//free space of vector
  return LIST;

}

void buildSets::freqofComslist( bool addOffset, int Koffset ){

  int i,j,K;

  Mmin=Clist[0].K;
  Mmax=Clist[0].K;

  if( addOffset ){//consider K starts from 0, add Koffset, i.e. 1.
    for( i=0; i<Clines; i++ ){ Clist[i].K = Clist[i].K + Koffset; }
  }
  
  for( i=0; i<Clines; i++ ){

    if( Clist[i].K > Mmax ) Mmax = Clist[i].K;
    if( Clist[i].K < Mmin ) Mmin = Clist[i].K;
    
  }
  
  //---unique list of communities and sizes
  for( i=Mmin; i<=Mmax; i++ ){
    COMS.push_back(pairIntInt(i,0));
  }

  for( i=0; i<Clines; i++ ){
    for( j=0; j<COMS.size(); j++ ){
      if( COMS[j].first == Clist[i].K){
	COMS[j].second++; continue;
      }	
    }
  }
 
 
  //Store number of communities
  M = COMS.size();  

}


LISTst* buildSets::freqofAnnolist( LISTst *_ALIST, int _ALINES, int &_F ){

  int i,j;

  vector<pairStrInt> ANNOIDS; 

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      ANNOIDS.push_back(pairStrInt(string(_ALIST[i].annoID),0));
    }
  }

  //---unique list of anno IDs & names
  sort(ANNOIDS.begin(),ANNOIDS.end(),sortpairStrInt());
  ANNOIDS.erase(unique(ANNOIDS.begin(),ANNOIDS.end()),ANNOIDS.end());

  _F = ANNOIDS.size();  
  LISTst *LIST = (LISTst*)calloc(_F,sizeof(LISTst));

  for( i=0; i<ANNOIDS.size(); i++ ){
    LIST[i].annoID = (char*)malloc((strlen(ANNOIDS[i].first.c_str())+1)*sizeof(char));
    strcpy(LIST[i].annoID, ANNOIDS[i].first.c_str());
    LIST[i].K  = 0;
    LIST[i].ID = 0;
  }

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      for( j=0; j<ANNOIDS.size(); j++ ){
	if( ANNOIDS[j].first.compare(_ALIST[i].annoID) == 0 ){
	  LIST[j].annoDES = (char*)malloc((strlen(_ALIST[i].annoDES)+1)*sizeof(char));
	  strcpy(LIST[j].annoDES, _ALIST[i].annoDES);
	  //LIST[j].K++;
	  //ANNOIDS[j].second++; continue;
	  continue;
	}	
      }

    }
  }

  vector<pairStrInt>().swap(ANNOIDS);  //free space of vector
  
  return LIST;
  

}


LISTst* buildSets::removeDuplicateIDs( LISTst *_ALIST, int &_ALINES, LISTst *_ANNOS, int _F ){

  int i,j;

  vector<pairStrInt> ANNOIDS; 

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      ANNOIDS.push_back(pairStrInt(string(_ALIST[i].annoID),_ALIST[i].ID));
    }
  }

  //---unique list of anno IDs & gene IDs
  sort(ANNOIDS.begin(),ANNOIDS.end(),sortpairStrInt());
  ANNOIDS.erase(unique(ANNOIDS.begin(),ANNOIDS.end()),ANNOIDS.end());
  
  //---delete _ALIST
  for(i=0; i<_ALINES; i++){
    if(_ALIST[i].annoID  != NULL) free(_ALIST[i].annoID);
    if(_ALIST[i].annoDES != NULL) free(_ALIST[i].annoDES);
  }
  free(_ALIST);
  _ALINES = 0;
  //---
  
  //---create new _ALIST struct
  _ALINES          = ANNOIDS.size();
  LISTst *ALISTtmp = (LISTst*)calloc(_ALINES,sizeof(LISTst));
  //---
  
  //loop through the set of annotation types, and fill _ALIST
  for( i=0; i<_ALINES; i++ ){

    for( j=0; j<_F; j++ ){

      if( ANNOIDS[i].first.compare(_ANNOS[j].annoID) == 0 ){

     ALISTtmp[i].annoID = (char*)malloc((strlen(_ANNOS[j].annoID)+1)*sizeof(char));
     strcpy(ALISTtmp[i].annoID, _ANNOS[j].annoID);
	
     ALISTtmp[i].annoDES = (char*)malloc((strlen(_ANNOS[j].annoDES)+1)*sizeof(char));
     strcpy(ALISTtmp[i].annoDES, _ANNOS[j].annoDES);
	
     ALISTtmp[i].ID = ANNOIDS[i].second; 	
     _ANNOS[j].K++;
	
      }

    }
  }

  vector<pairStrInt>().swap(ANNOIDS);  //free space of vector

  return(ALISTtmp);
  
}
  
buildSets::~buildSets(){ 

  int i,j;

  freeSpace();

  if( Clist != 0 && freedClist == false ){ 
    for(i=0; i<Clines; i++){ 
      if(Clist[i].annoID != NULL) free(Clist[i].annoID);
      if(Clist[i].annoDES != NULL)free(Clist[i].annoDES);      
    }
    free(Clist);
    freedClist = true;
  }

    
}

void buildSets::assignSpace(){}

void buildSets::freeSpace(){

 int i,j;

  if( freedfiles == false ){
    for( i=0; i<files.size(); i++){ 
      delete files[i]; vector<fileReader*>().swap(files); freedfiles=true; }
  }   

  if( freedAlist == false ){
    for(i=0; i<Alines.size(); i++){
      for(j=0; j<Alines[i]; j++){
	if(Alist[i][j].annoID  != NULL) free(Alist[i][j].annoID);
	if(Alist[i][j].annoDES != NULL) free(Alist[i][j].annoDES);
      }
      free(Alist[i]);
    }
    freedAlist = true;
  }

  if( freedANNOS == false ){
    for(i=0; i<Fsize.size(); i++){
      for(j=0; j<Fsize[i]; j++){
	if(ANNOS[i][j].annoID  != NULL) free(ANNOS[i][j].annoID);
	if(ANNOS[i][j].annoDES != NULL) free(ANNOS[i][j].annoDES);
      }
      free(ANNOS[i]);
    }
    freedANNOS = true;
  }
  
  
}

void buildSets::printList( LISTst *LIST, int SIZE){
  
  int i;
  
  for(i=0; i<SIZE; i++){
    printf("%s\t",LIST[i].annoID); 
    printf("%s\t",LIST[i].annoDES); 
    cout << LIST[i].ID      << "\t"
         << LIST[i].K       << endl;
    
  }

}

int buildSets::getKMin(){ return Mmin; }

void buildSets::addKOffset( int addOffset ){

  M_old    = M;
  Mmin_old = Mmin;
  Mmax_old = Mmax;

  COMS.clear();
  
  if( addOffset < 0 ){ addOffset = 0; }
  
  //consider K = 0
  freqofComslist( true, addOffset );

  /*
  cout << "------------------- "<< endl;
  cout << "Before resetting"    << endl;
  cout << "M    : " << M_old    << endl;
  cout << "K min: " << Mmin_old << endl;
  cout << "K max: " << Mmax_old << endl;

  cout << "After resetting" << endl;
  cout << "M    : " << M    << endl;
  cout << "K min: " << Mmin << endl;
  cout << "K max: " << Mmax << endl;
  cout << "------------------- "<< endl;
  */
  
}

