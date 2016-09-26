#include "pc.h"
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
#include <zlib.h>

extern "C" {
// pliu 20160911
//#include "R.h"
//#include "Rmath.h"
//////
#include "Rinternals.h"
#include "Rdefines.h"
#include "maqmap.h"
}

using namespace std;
using namespace __gnu_cxx; 


class lessAbsoluteValue {
public:
  bool operator()(int a, int b) const {
    return abs(a) < abs(b);
  }
};



//#define DEBUG 1

extern "C" {

  // read in text version of maq map
  SEXP read_binmaqmap(SEXP filename,SEXP read_tag_names_R) {

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
  

  gzFile f=gzopen(fname,"r");

  maqmap_t *m = maqmap_read_header(f);  
  maqmap1_t *m1, mm1;
  m1 = &mm1;

  if (!f)  { 
    cout<<"can't open input file \""<<fname<<"\"\n"; 
  }  else {
    Rprintf("opened %s\n",fname);

    // read in bed line
    string line;
    int fcount=0;
    while(maqmap_read1(f, m1)) {
      string tagname=string(m1->name);
      string chr=string(m->ref_name[m1->seqid]);
      int len=m1->size;
      int fpos=(m1->pos>>1) + 1;
      if(m1->pos&1) {
	fpos=-1*(fpos+len-1);
      }
      int nm=m1->info1&0xf;

#ifdef DEBUG  
      Rprintf("read in map line chr=%s tagname=%s fpos=%d, nm=%d, len=%d\n",chr.c_str(),tagname.c_str(),fpos,nm,len);
#endif
    

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
    gzclose(f);
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


}
