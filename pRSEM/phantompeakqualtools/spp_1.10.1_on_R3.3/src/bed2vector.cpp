#include "pc.h"
#include "config.h"
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <algorithm>
#include <string>
#include <functional>
#include <utility>
#include <ext/hash_map>
#include <boost/tokenizer.hpp>

#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif 

extern "C" {
// pliu 20160911
//#include "R.h"
//#include "Rmath.h"
//////
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 


class lessAbsoluteValue {
public:
  bool operator()(int a, int b) const {
    return abs(a) < abs(b);
  }
};



#ifdef HAVE_LIBBZ2
int get_bzline(BZFILE* b,string& line) {
  char c;
  int     nBuf;
  int bzerror=BZ_OK;

  while(bzerror == BZ_OK)  {  
    nBuf=BZ2_bzRead(&bzerror, b, &c, 1);
    if(bzerror==BZ_OK) {
      if(c=='\n') {
	return bzerror;
      } else {
	line+=c;
      }
    }
  }
  return bzerror;
}

int get_a_line(FILE *f,BZFILE *b,int bz2file,string& line) {
  line="";
  if(bz2file) {
    int bzerror=get_bzline(b,line);
    if(bzerror==BZ_OK) {
      return(1);
    } else {
      if(bzerror!=BZ_STREAM_END) {
	cerr<<"encountered BZERROR="<<bzerror<<endl;
      }
      return(0);
    }
  } else {
    char *cline=NULL;
    size_t n;
    if(getline(&cline,&n,f) != -1) {
      if(cline) {
	cline[strlen(cline)-1]='\0';
	line+=cline;
	free(cline);
      }
      return(1);
    } else {
      return(0);
    }
  }
}
#endif


/**
 * Read in .bed data into a list chromosome of vectors representing 5' positions, with sign
 * corresponding to the strand.
 */

//#define DEBUG 1

extern "C" {
SEXP read_bed_ends(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");


  ifstream bed_file(fname);

#ifdef DEBUG  
  Rprintf("opened %s\n",fname);
#endif

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
    
  int fcount=0;
  while(getline(bed_file,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++; //chr=chr.substr(3,strlen(chr.c_str()));
      string str_start=*sit++;
      int fstart=atoi(str_start.c_str());
      string str_end=*sit++;
      int fend=atoi(str_end.c_str());
      int fpos=fstart;
      if(sit!=tok.end()) {
         string u0=*sit++;
         string nfield=*sit++;
         string strand=*sit++;
         if(strand=="-") { 
	   fpos=-1*fend;
         }
      }

      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d\n",chr.c_str(),cind,fpos);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  bed_file.close();
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
    sort(csi->begin(), csi->end(), lessAbsoluteValue());
  }

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    SEXP nv;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_nv=INTEGER(nv);
    int i=0;
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_nv[i++]=*pi;
    }
    SET_VECTOR_ELT(ans, csi-pos.begin(), nv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



SEXP read_meland_old(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");


  ifstream bed_file(fname);

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
    
  int fcount=0;
  while(getline(bed_file,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      sit++; sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string str_len=*sit++;
      int len=atoi(str_len.c_str());
      string chr=*sit++; chr=chr.substr(3,strlen(chr.c_str()));
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  bed_file.close();
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    
    
    
    SEXP tv,nv,lv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 3));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  int get_a_line(FILE *f,string& line) {
    line="";
    char cline[1024];
    if(fgets(cline,1024,f)) {
      line+=cline;
      return(1);
    } else {
      return(0);
    }
  }


  SEXP read_meland(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  
  Rprintf("opened %s\n",fname);


  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string str_len=*sit++;
      int len=atoi(str_len.c_str());
      string chr=*sit++; chr=chr.substr(3,strlen(chr.c_str()));
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 3, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,lv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 3+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 3, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



// reads regular eland files, recording mismatch positions
SEXP read_eland_mismatches(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > mm1; // position of the first mismatch (or 0 for none)
  vector< vector<int> > mm2; // position of the second mismatch

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      sit++; 
      string seq=*sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string chr=*sit++; 
      // extract chromosome name from this
      int chrp=chr.find("chr");
      int pp=chr.find('.');
      chr=chr.substr(chrp+3,pp-chrp-3);
      
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());


      string strand=*sit++;
      int nstrand=0;
      if(strand=="R") { 
	fpos=-1*(fpos+seq.size()-1);
	nstrand=1;
      }

      sit++;
      
      int nm1=0; int nm2=0;
      if(sit!=tok.end()) {
	string nms=*sit++;
	nm1=atoi(nms.substr(0,nms.size()-1).c_str());
	if(nstrand) { nm1=seq.size()-nm1+1; }
      }
      if(sit!=tok.end()) {
	string nms=*sit++;
	nm2=atoi(nms.substr(0,nms.size()-1).c_str());
	if(nstrand) { nm2=seq.size()-nm2+1; }
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	mm1.push_back(vector<int>());
	mm2.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (mm1[cind]).push_back(nm1);
      (mm2[cind]).push_back(nm2);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm1=%d, nm2=%d\n",chr.c_str(),cind,fpos,nm1,nm2);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=mm1.begin()+(csi-pos.begin());
    lsi=mm2.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("f"));
    SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    
    
    
    SEXP tv,nv,lv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 3));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in regular eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_eland(SEXP filename,SEXP read_tag_names_R,SEXP eland_tag_length_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
  int eland_tag_length=*(INTEGER(eland_tag_length_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string sequence=*sit++;
      int len=sequence.size();
      // adjust probe length if eland length limit was specified
      if(eland_tag_length>0 && len>eland_tag_length) {
	len=eland_tag_length;
      }
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string chr=*sit++; 
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;

      if(str_strand[0]=='R') {
	fpos=-1*(fpos+len-1);
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



  // read in extended eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_eland_extended(SEXP filename,SEXP read_tag_names_R,SEXP eland_tag_length_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
  int eland_tag_length=*(INTEGER(eland_tag_length_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string machinename=*sit++;
      string runnumber=*sit++;
      string lanenumber=*sit++;
      *sit++;
      
      string str_x=*sit++;
      string str_y=*sit++;

      string tagname=machinename+"."+runnumber+"."+lanenumber+"."+str_x+"."+str_y;

      

      *sit++;
      *sit++;

      
      string sequence=*sit++;
      *sit++;
      
      string chr=*sit++; 
      string contig=*sit++; 
      chr=chr+contig;
      
      int len=sequence.size();
      // adjust probe length if eland length limit was specified
      if(eland_tag_length>0 && len>eland_tag_length) {
	len=eland_tag_length;
      }


      
      string str_pos=*sit++;
      if(str_pos.size()<1) { continue; }
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;

      if(str_strand[0]=='R') {
	fpos=-1*(fpos+len-1);
      }

      string str_nm=*sit++;
      // count non-digit characters
      int nm=0;
      for(int i=0;i<str_nm.size();i++) {
	if(!isdigit(str_nm[i])) { nm++; }
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in eland multi files, adjusting the negative strand coordinate by sequence length
SEXP read_eland_multi(SEXP filename,SEXP read_tag_names_R,SEXP eland_tag_length_R) {
  
#ifdef DEBUG  
  Rprintf("read_eland_muti() : start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
  int eland_tag_length=*(INTEGER(eland_tag_length_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t","");
  boost::char_separator<char> comsep(",","",boost::keep_empty_tokens);
  boost::char_separator<char> colsep(":","",boost::keep_empty_tokens);
  
  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int nline=0;
  int fcount=0;
  while(get_a_line(f,line)) {
    nline++;
    // chomp
    size_t elpos = line.find_last_not_of("\n");
    if(elpos != string::npos) {
      line = line.substr(0, elpos+1);
    }
#ifdef DEBUG  
    Rprintf("line %d: %s\n",nline,line.c_str());
#endif

    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string sequence=*sit++;
      string mspec=*sit++;
      // parse out match spec
      
      if(mspec=="NM" || mspec=="QC") { continue; }
#ifdef DEBUG  
      Rprintf("parsing out spec \"%s\" : ",mspec.c_str());
#endif
      
      tokType stok(mspec, colsep);
      tokType::iterator ssit=stok.begin();
      string str_nm0=*ssit++;
      
      int nm=0;
      int nm0=atoi(str_nm0.c_str());
      if(nm0>1) { 
#ifdef DEBUG  
	Rprintf("rejected for nm0\n");
#endif
	continue; 
      }
      if(nm0==0) {
	string str_nm1=*ssit++;
	int nm1=atoi(str_nm1.c_str());
	if(nm1>1) { 
#ifdef DEBUG  
	  Rprintf("rejected for nm1\n");
#endif
	  continue; 
	}
	if(nm1==0) {
	  string str_nm2=*ssit++;
	  int nm2=atoi(str_nm2.c_str());
	  if(nm2>1) { 
#ifdef DEBUG  
	    Rprintf("rejected for nm2\n");
#endif
	    continue; 
	  }
	  nm=2;
	} else {
	  nm=1;
	}
      }

#ifdef DEBUG  
      Rprintf("accepted (nm=%d)\n",nm);
#endif
      int npos=0;
      string mpos=*sit++;
      vector<string> mposc;
      vector<int> mposp;
      tokType ptok(mpos, comsep);
      string prevchr;
      for(tokType::iterator psit=ptok.begin();psit!=ptok.end();psit++) {
	string cpos=*psit;
	npos++;
	int strand=1;
	if(cpos.size()<5) {
	  Rprintf("ERROR: line=%d, match %d is too short: \"%s\"; ",nline,npos,cpos.c_str());
	}
	char lc=cpos.at(cpos.size()-1);
	
	if(atoi(&lc)==nm) {
	  switch(cpos.at(cpos.size()-2)) {
	  case 'R': strand=-1; break;
	  case 'F': strand=1; break;
	  default:
	    Rprintf("ERROR: line=%d, match %d specifies an invalid strand %c\n",nline,npos,cpos.at(cpos.size()-2)); break;
	    continue;
	  }
          string chr,str_pos;
	  size_t colpos=cpos.find(":");
	  if(colpos==string::npos) {
            if(npos>1) {
              chr=prevchr;
              str_pos=cpos.substr(0,cpos.size()-2);
            } else {
	      Rprintf("ERROR: line=%d, match %d does not contain chromosome separator: \"%s\"\n",nline,npos,cpos.c_str()); 
	      continue;
            }
	  } else {
	      chr=cpos.substr(0,colpos);
	      str_pos=cpos.substr(colpos+1,cpos.size()-3-colpos);
          }
#ifdef DEBUG  
	  Rprintf("\"%s\" : chr=%s, pos=%s, strand=%d\n",cpos.c_str(),chr.c_str(),str_pos.c_str(),strand);
#endif	  
	  int pos=strand*atoi(str_pos.c_str());
	  mposc.push_back(chr);
	  mposp.push_back(pos);
	}
      }

      string chr;
      int fpos;
      if(mposc.size()!=1) {
	if(mposc.size()==0) {
	  Rprintf("ERROR: line=%d: no %d-mismatch matches were found in \"%s\"\n",nline,nm,mpos.c_str()); 
	} else {
	  Rprintf("ERROR: line=%d: more than one (%d) %d-mismatch matches were found in \"%s\"\n",nline,mposc.size(),nm,mpos.c_str()); 
	}
	continue;
      } else {
	chr=*mposc.begin();
	fpos=*mposp.begin();
      }
      
      int len=sequence.size();
      // adjust probe length if eland length limit was specified
      if(eland_tag_length>0 && len>eland_tag_length) {
	len=eland_tag_length;
      }

      if(fpos<0) {
	fpos=-1*(-1*fpos+len-1);
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in regular eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_bowtie(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);
  boost::char_separator<char> sep2(",");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; 
  } else {
#ifdef HAVE_LIBBZ2
    BZFILE* b;  
    int bzerror;
    
    int bz2file=0;
    if(strstr(fname,".bz2")) {
      bz2file=1;
      b=BZ2_bzReadOpen (&bzerror, f, 0, 0, NULL, 0);
      if (bzerror != BZ_OK)  { cout<<"bzerror="<<bzerror<<endl; }
    }
#endif

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
#ifdef HAVE_LIBBZ2
  while(get_a_line(f,b,bz2file,line)) {
#else
  while(get_a_line(f,line)) {
#endif

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string str_strand=*sit++;
      string chr=*sit++; 

      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());

      string sequence=*sit++;
      sit++; sit++;
      string mm=*sit++;

      int len=sequence.size();
      if(str_strand[0]=='-') {
	fpos=-1*(fpos+len-1);
      }
      // determine number of mismatches
      int nm=0;
      if(mm.size()>0) {
	nm++;
	string::size_type tp(0);
	while(tp!=string::npos) {
	  tp = mm.find(",",tp);
	  if(tp!=string::npos) {
	    tp++;
	    ++nm;
	  }
	}
      }


      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }

#ifdef HAVE_LIBBZ2
  BZ2_bzReadClose( &bzerror, b);
#endif
 fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in helicos tab-separated alignment output (regular or bz2)
  SEXP read_helicostabf(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length of the match
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);
  boost::char_separator<char> sep2(",");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; 
  } else {
#ifdef HAVE_LIBBZ2
    BZFILE* b;  
    int bzerror;
    
    int bz2file=0;
    if(strstr(fname,".bz2")) {
      bz2file=1;
      b=BZ2_bzReadOpen (&bzerror, f, 0, 0, NULL, 0);
      if (bzerror != BZ_OK)  { cout<<"bzerror="<<bzerror<<endl; }
    }
#endif

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  int nlines=0;
#ifdef HAVE_LIBBZ2
  while(get_a_line(f,b,bz2file,line)) {
#else
  while(get_a_line(f,line)) {
#endif

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif
    nlines++;
    // skip comments
    if(line[0]=='#') { continue; }
    if(line.compare(0,12,"Reference_ID")==0) { 
#ifdef DEBUG  
      Rprintf("matched header on line %d\n",nlines); 
#endif
      continue; 
    }

    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++; 
      string tagname=*sit++;
      string str_startpos=*sit++;
      string str_endpos=*sit++;

      string str_tstart=*sit++;
      string str_tend=*sit++;
      int len=atoi(str_tend.c_str())-atoi(str_tstart.c_str());

      sit++; sit++;
      string str_ndel=*sit++;
      string str_nins=*sit++;
      string str_nsub=*sit++;
      
      string str_strand=*sit++;
      int fpos;
      if(str_strand[0]=='-') {
	fpos=-1*atoi(str_endpos.c_str()); 
      } else {
	fpos=atoi(str_startpos.c_str()); 
      }

      // determine number of mismatches
      int nm=atoi(str_ndel.c_str())+atoi(str_nins.c_str())+atoi(str_nsub.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d\n",chr.c_str(),cind,fpos,nm);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }

#ifdef HAVE_LIBBZ2
  BZ2_bzReadClose( &bzerror, b);
#endif
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<int> >::const_iterator lsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 3, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,lv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator lni=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*lni++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 3+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 3, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



  // read in text version of maq map
  SEXP read_maqmap(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string chr=*sit++;
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;
      sit++; sit++; sit++; sit++; sit++; 
      string str_nm=*sit++;
      sit++; sit++; sit++; 
      string str_len=*sit++;
      int nm=atoi(str_nm.c_str());
      int len=atoi(str_len.c_str());

      if(str_strand[0]=='-') {
	fpos=-1*(fpos+len-1);
      }

      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}





  // read in tagalign file
  SEXP read_tagalign(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++;
      string str_spos=*sit++;
      string str_epos=*sit++;
      sit++; 
      string str_qual=*sit++;
      string str_strand=*sit;

      int fpos;
      if(str_strand[0]=='+') {
	fpos=atoi(str_spos.c_str());
      } else {
	fpos=-1*atoi(str_epos.c_str());
      }
      int nm=atoi(str_qual.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d nm=%d\n",chr.c_str(),cind,fpos,nm);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    
    
    SEXP tv,nv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 2));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}




  // arachne madness
  SEXP read_arachne(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  



  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {

#ifdef HAVE_LIBBZ2
    BZFILE* b;  
    int bzerror;
    
    int bz2file=0;
    if(strstr(fname,".bz2")) {
      bz2file=1;
      b=BZ2_bzReadOpen (&bzerror, f, 0, 0, NULL, 0);
      if (bzerror != BZ_OK)  { cout<<"bzerror="<<bzerror<<endl; }
    }
#endif


  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
#ifdef HAVE_LIBBZ2
  while(get_a_line(f,b,bz2file,line)) {
#else
  while(get_a_line(f,line)) {
#endif

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++;
      string str_spos=*sit++;
      int nm=0;
      if(sit!=tok.end()) {
	string str_mm=*sit;
	nm=atoi(str_mm.c_str());
      }
      
      int fpos=atoi(str_spos.c_str());;
      
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d nm=%d\n",chr.c_str(),cind,fpos,nm);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
#ifdef HAVE_LIBBZ2
  BZ2_bzReadClose( &bzerror, b);
#endif

  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    
    
    SEXP tv,nv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 2));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // arachne madness
  SEXP read_arachne_long(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length of the match

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  



  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {

#ifdef HAVE_LIBBZ2
    BZFILE* b;  
    int bzerror;
    
    int bz2file=0;
    if(strstr(fname,".bz2")) {
      bz2file=1;
      b=BZ2_bzReadOpen (&bzerror, f, 0, 0, NULL, 0);
      if (bzerror != BZ_OK)  { cout<<"bzerror="<<bzerror<<endl; }
    }
#endif


  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
#ifdef HAVE_LIBBZ2
  while(get_a_line(f,b,bz2file,line)) {
#else
  while(get_a_line(f,line)) {
#endif

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string query=*sit++;
      if(query!="QUERY") { continue; }
      *sit++; *sit++; *sit++; *sit++; 
      string str_strand=*sit++;
      string chr=*sit++;
      string str_startpos=*sit++;
      string str_endpos=*sit++;
      
      int fpos;
      if(str_strand[0]=='1') {
	fpos=-1*atoi(str_endpos.c_str()); 
      } else {
	fpos=atoi(str_startpos.c_str()); 
      }
#ifdef DEBUG  
      Rprintf("chr=%s, fpos=%d\n",chr.c_str(),fpos);
#endif
      *sit++;
      string str_nblocks=*sit++;
      int nblocks=atoi(str_nblocks.c_str());
#ifdef DEBUG  
      Rprintf("nblocks=%d\n",nblocks);
#endif
      // tally up the read length and the number of mismatches for all blocks
      int len=0; int nm=0;
      for(int i=0;i<nblocks;i++) {
	string str_sgs=*sit++;
	int sgs=atoi(str_sgs.c_str());
	string str_slen=*sit++;
	int slen=atoi(str_slen.c_str());
	string str_snm=*sit++;
	int snm=atoi(str_snm.c_str());
#ifdef DEBUG  
	Rprintf("sgs=%d, slen=%d, snm=%d\n",sgs,slen,snm);
#endif
	len+=slen;
	nm+=abs(sgs)+snm;
      }
      nm+=nblocks-1;
      
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d nm=%d len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
#ifdef HAVE_LIBBZ2
  BZ2_bzReadClose( &bzerror, b);
#endif

  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<int> >::const_iterator lsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    
    
    SEXP tv,nv,lv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator lni=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*lni++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 3));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


}
