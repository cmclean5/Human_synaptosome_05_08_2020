#include "Headers.h"

//define guards, so headers are declare only once.
#ifndef FILEREADER_H
#define FILEREADER_H

class fileReader {

 public:
  fileReader();
  fileReader(FILE *);
 ~fileReader();
 
 int  read_line   (char *);
 int  read_line   (char *, char);
 int  fill_buffer ();
 void free_buffer ();
 void reset_buffer();
 int  next_line   (char *);
 int  get_nlines  ();
 char*** readFile ( int );


static inline void subString(char s[], char sub[], int p, int l){
   int c    = 0;
   int diff = (l-p); 
 
   while (c <= diff) {
     sub[c] = s[p+c-1];
     c++;
   }
   sub[c] = '\0';
 };
 

static inline void findDelimiterPos( char Line[], int Length, char dels[], int Ndels, vector<int> &pos){ 

  int i,d;
  
  //loop over each character in Line
  for(i=0; i<Length; i++){
      
    //search position of delemiters in line
    for(d=0; d<Ndels; d++){
      if(Line[i]==dels[d]){ pos.push_back(i); }
    }

  }

};


static inline int readDirectory( const char *Dir, vector<string> &FILES ){

  string filepath;
  DIR *dp;
  struct dirent *dirp;
  struct stat filestat;

  dp = opendir( Dir );
  if (dp == NULL)
    {
      cout << "Error(" << errno << ") opening " << Dir << endl;
      return errno;
    }

  while ((dirp = readdir( dp )))
    {
      char* current_file = dirp->d_name;
      
      filepath = string(Dir) + "/" + current_file;
      
      // If the file is a directory (or is in some way invalid) we'll skip it 
      if (stat( filepath.c_str(), &filestat )) continue;
      if (S_ISDIR( filestat.st_mode ))         continue;
      
      if( current_file[0]=='.' )               continue;
      
      // save file into list
      FILES.push_back( filepath );
      
    }
  
  closedir( dp );
  
  return 0;

};


 static inline int writeDirectory( const char *Dir ){

   int status;
   char buffer[250];

   if( Dir==NULL ) return -1;

   status = mkdir(Dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

   return status;

 };

//Get current date/time
static inline void currentDateTime( string &result) {
    time_t now = time(0);
    char *dt   = ctime(&now);
    int LENGTH = strlen(dt);
    dt[LENGTH-1] = '\0';

    result     = string(dt)+"/";
    
};

 private:
 
 int Nlines;
 int Ncols;

 FILE *stream;

 char ***DATASET;

 // Constants
 #define LINELENGTH 1000

 // Types
 struct LINE {
   char *str;
   struct LINE *ptr;
 };

 // Globals
 LINE *first;
 LINE *current;

 //set file delineator(s)
  static const int DELSIZE = 1;
  char dels[DELSIZE];

  //set file header delineator(s)
  static const int HEADDELSIZE = 1;
  char HEADdel[HEADDELSIZE];
  
};

#endif
