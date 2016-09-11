#include <math.h>
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"


#undef DEBUG 1

// dout is npos-length output array.
// n - number of positions in pos (and length of tc count array)
// spos - starting position
void cdensum(int *n, double *pos, double *tc, double *spos, int *bw,int *dw, int *npos, int *step,double *dout)
{
  int i,j;
 
  double epos= *spos + ((double) *npos);
  double dbw=(double) *bw;
  for(i = 0; i< *n; i++) {
    // size of the window to which the contributions should be added
    int in=(int) (pos[i]- *spos);
    int ic=tc[i];
    int whs=(*dw)*(*bw)*ic;
    int ws=(int) floor((in-whs)/(*step));
    int we=(int) ceil((in+whs)/(*step));
    if(ws<0) { ws=0; } 
    if(we>= *npos) { we= *npos -1; }
    
    for(j=ws;j<we;j++) {
      double beta=((double)(j*(*step)-in))/dbw;
      dout[j]+=((double)ic)*exp(-0.5*beta*beta);
    }
  }
}


// window tag counts
// dout is npos-length output array that will contain window tag counts
// windows are of a specified size, moved at a specified step
// n - number of positions in sorted tag array (positive only)
// spos - starting position
void window_n_tags(int *n, double *pos, double *spos, int *window_size, int *window_step, int *npos, int *dout)
{
  int i;
  int cs=0; int ce=0; // current array start/end indecies
  int ctc=0; // current tag count
  double wpos=*spos-(*window_size)/2; // left-edge position
  //Rprintf("n=%d; window_size=%d, window_step=%d, npos=%d, spos=%f\n",*n,*window_size,*window_step,*npos,*spos);
  for(i=0;i<*npos;i++) {
    // advance end if needed
    double ep=wpos+(*window_size);
    while(ce<(*n) && pos[ce]<=ep) {
      ctc++; ce++;
    }
    // advance start
    while(cs<*n && pos[cs]<wpos) {
      ctc--; cs++;
    }
    dout[i]=ctc;
    // advance window position
    wpos+=*window_step;
  }
}

// window tag counts
// windows are of a specified size, moved at a specified step
// pos - tag positions (positive, pre-shifted)y
// spos - starting position
// returns nsteps-length output array that will contain window tag counts
SEXP cwindow_n_tags(SEXP pos_R, SEXP spos_R, SEXP window_size_R, SEXP window_step_R, SEXP nsteps_R) {
  double* pos=REAL(pos_R);
  int n=LENGTH(pos_R);
  int window_size=*INTEGER(window_size_R);
  int window_step=*INTEGER(window_step_R);
  int nsteps=*INTEGER(nsteps_R);
  double spos=*REAL(spos_R);
  
  // allocate return array
  SEXP tc_R;
  PROTECT(tc_R=allocVector(INTSXP,nsteps));
  int* dout=INTEGER(tc_R);

  int i;
  int cs=0; int ce=0; // current array start/end indecies
  int ctc=0; // current tag count
  double wpos=spos-window_size/2; // left-edge position
  //Rprintf("n=%d; window_size=%d, window_step=%d, npos=%d, spos=%f\n",n,window_size,window_step,nsteps,spos);
  for(i=0;i<nsteps;i++) {
    // advance end if needed
    double ep=wpos+window_size;
    while(ce<n && pos[ce]<=ep) {
      ctc++; ce++;
    }
    // advance start
    while(cs<n && pos[cs]<wpos) {
      ctc--; cs++;
    }
    dout[i]=ctc;
    // advance window position
    wpos+=window_step;
  }
  UNPROTECT(1);
  return(tc_R);
}

// tag counts in windows around specified positions
// pos - tag positions 
// ntags - number of tags in each position
// wpos - window positions
// returns a pos-length vector giving number of tags that fall within window_half_size from the provided positions
SEXP cwindow_n_tags_around(SEXP pos_R, SEXP ntags_R, SEXP wpos_R, SEXP window_half_size_R) {
  double* pos=REAL(pos_R);
  int* ntags=INTEGER(ntags_R);
  int n=LENGTH(pos_R);
  double* wpos=REAL(wpos_R);
  int nw=LENGTH(wpos_R); // number of windows
  double whs=(double) *INTEGER(window_half_size_R);
  
  // allocate return array
  SEXP tc_R;
  PROTECT(tc_R=allocVector(INTSXP,nw));
  int* dout=INTEGER(tc_R);

  int i;
  int cs=0; int ce=0; // current array start/end indecies
  int ctc=0; // current tag count
  for(i=0;i<nw;i++) {
    //if(i>(nw-2)) {      Rprintf("-i=%d; cs=%d, ce=%d; ctc=%d\n",i,cs,ce,ctc);    }    
    // advance end if needed
    double ep=wpos[i]+whs;
    while(ce<n && pos[ce]<=ep) {
      ctc+=ntags[ce]; ce++;
    }
    // advance start
    double sp=wpos[i]-whs;
    while(cs<n && pos[cs]<sp) {
      ctc-=ntags[cs]; cs++;
    }
    dout[i]=ctc;
    // if(i>(nw-2)) {      Rprintf("+i=%d; cs=%d, ce=%d; ctc=%d\n",i,cs,ce,ctc);    }
  }
  UNPROTECT(1);
  return(tc_R);
}

