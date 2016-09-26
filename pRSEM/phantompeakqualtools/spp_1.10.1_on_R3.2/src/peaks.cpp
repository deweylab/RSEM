#include <vector>
#include <string.h>
#include <iostream>
#include <string>
#include <set>

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 

/**
 * Calculate all local peaks
 */

//#define DEBUG 1

extern "C" {
  SEXP find_peaks(SEXP x_R,SEXP thr_R,SEXP max_span_R) {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    double* x=REAL(x_R);
    int nx=LENGTH(x_R);
    int max_span=*INTEGER(max_span_R);
    double thr=REAL(thr_R)[0];
#ifdef DEBUG  
    Rprintf("n=%d; thr=%f; max_span=%d\n",nx,thr,max_span);
#endif

    vector<int> pos;
  
    double pv=x[0];
    double ppv=0; // previous peak value
    int ppp=-max_span-1; // previous peak position
    
    for(int i=1;i<(nx-1);i++) {
      if(x[i]>pv && x[i]>=thr && x[i]>x[i+1]) {
	if(max_span>2) {
	  //Rprintf("i=%d; ppp=%d\n",i,ppp);
	  if(i-ppp > max_span) {
	    if(ppp>=0) {
	      pos.push_back(ppp);
	    }
	    //Rprintf("recorded %d; now %d\n",ppp,i);
	    ppp=i; ppv=x[i];
	  } else {
	    if(x[i]>ppv) {
	      //Rprintf("reset from %d to %d\n",ppp,i);
	      ppp=i; ppv=x[i];
	    }
	  }
	} else {
	  pos.push_back(i);
	}
      }
      if(x[i]!=x[i+1]) { pv=x[i]; }
    }

    // add remaining peak
    if(max_span>2 && ppp>=0) {
      pos.push_back(ppp);
    }

    SEXP nv;
    PROTECT(nv=allocVector(INTSXP,pos.size())); 
    int* i_nv=INTEGER(nv);
    int i=0;
    for(vector<int> ::const_iterator pi=pos.begin();pi!=pos.end();++pi) {
      i_nv[i++]=1+(*pi);
    }
  
    UNPROTECT(1);
    return(nv);
  }




  /************************************************************************/
  // given a data vector d (positive values) and a set of signed center coordinates pos,
  // returns coordinates of data points relative to the centers
  // size is the size of the region around the centers
  // return: vector of relative coordinates (x) and indecies of centers relative the coordinate
  // was calculated (i).
  SEXP get_relative_coordinates(SEXP d_R,
				SEXP pos_R,
				SEXP size_R)
  {
    int *d, *pos; 
    int npos,nd,size;
  
    d = INTEGER(d_R); pos = INTEGER(pos_R);
    npos=LENGTH(pos_R);  nd=LENGTH(d_R);
    size = INTEGER(size_R)[0];
#ifdef DEBUG  
    Rprintf("|d|=%d, |c|=%d, size=%d\n",nd,npos,size);
#endif

    vector<int> x; vector<int> xi;
    int k=0; // current pos index
    
    for(int i=0;i<nd;i++) {
      // increment k until pos[k]+size>=d[i]
      while((abs(pos[k])+size) < d[i]) { k++; if(k==npos) { break; };
#ifdef DEBUG  
	Rprintf("advancing k to %d\n",k);
#endif
      }
      if(k==npos) { break; };
      // increment i until d[i]>=pos[k]-size
      while((abs(pos[k])-size) > d[i]) { i++; if(i==nd) { break; }
#ifdef DEBUG  
	Rprintf("advancing i to %d\n",i);
#endif
      }
      if(i==nd) { break; }


      int l=k;
      while((l<npos) && ((abs(pos[l])-size) <= d[i])) { l++; 
#ifdef DEBUG  
	Rprintf("advancing l to %d\n",l);
#endif
      }
      for(int j=k;j<l;j++) {
	int pd=d[i]-abs(pos[j]);
	if(abs(pd)<=size) {
	  // record
	  if(pos[j]>0) {
	    x.push_back(pd);
	  } else {
	    x.push_back(-1*pd);
	  }
	  xi.push_back(j);
#ifdef DEBUG  	
	  Rprintf("recorded i=%d, j=%d\n",i,j);
#endif
	} else {
	  break;
	}
      }
    }
    
    SEXP xv_R,xiv_R;
    PROTECT(xv_R=allocVector(INTSXP,x.size())); 
    PROTECT(xiv_R=allocVector(INTSXP,x.size())); 
    int* xv=INTEGER(xv_R);
    int* xiv=INTEGER(xiv_R);

    int i=0;
    for(vector<int> ::const_iterator pi=x.begin();pi!=x.end();++pi) {
      xv[i++]=*pi;
    }
    i=0;
    for(vector<int> ::const_iterator pi=xi.begin();pi!=xi.end();++pi) {
      xiv[i++]=1+(*pi);
    }
    
    SEXP ans_R, names_R;
    PROTECT(names_R = allocVector(STRSXP, 2));
    SET_STRING_ELT(names_R, 0, mkChar("x"));
    SET_STRING_ELT(names_R, 1, mkChar("i"));
    
    PROTECT(ans_R = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans_R, 0, xv_R);
    SET_VECTOR_ELT(ans_R, 1, xiv_R);
    setAttrib(ans_R, R_NamesSymbol, names_R);
  
    UNPROTECT(4);
    return(ans_R);
  }


  // determines a set of points within a set of fragments
  // note: all vectors sorted in ascending order
  // note: all vectors are integers
  // x_R - vector of point positions 
  // se_R - vector of start and end positions 
  // fi_R - vector of signed fragment indecies
  // return_list_R - whether a list of fragments should be returned for each point
  // return_unique_R - whether points in multiple fragments should be omitted
  SEXP points_within(SEXP x_R,SEXP se_R,SEXP fi_R,SEXP return_list_R,SEXP return_unique_R,SEXP return_point_counts_R) {
#ifdef DEBUG  
    Rprintf("start\n");
#endif
    int* x=INTEGER(x_R);
    int nx=LENGTH(x_R);
    int* se=INTEGER(se_R);
    int* fi=INTEGER(fi_R);
    int nf=LENGTH(se_R);

    int return_list=*(INTEGER(return_list_R));
    int return_unique=*(INTEGER(return_unique_R));
    int return_point_counts=*(INTEGER(return_point_counts_R));

#ifdef DEBUG  
    Rprintf("nf=%d; nx=%d, return_list=%d, return_unique=%d, return_point_counts=%d\n",nf/2,nx,return_list,return_unique,return_point_counts);
#endif
    set<int> fset;


    SEXP nv; int *i_nv;
    int np=0;
    if(return_point_counts) {
      PROTECT(nv = allocVector(INTSXP, nf/2)); np++;      
      i_nv=INTEGER(nv);
      for(int i=0;i<nf/2;i++) { i_nv[i]=0; }
    } else if(return_list) {
      PROTECT(nv = allocVector(VECSXP, nx)); np++;
    } else {
      PROTECT(nv=allocVector(INTSXP,nx));  np++;
      i_nv=INTEGER(nv);
    }

    int j=0;

    for(int i=0;i<nx;i++) {
      // advance j
      while(j<nf && se[j]<x[i]) {
	int frag=fi[j];
	if(frag>0) { // insert
	  fset.insert(frag);
#ifdef DEBUG  
	  Rprintf("inserted frag %d, size=%d\n",frag,fset.size());
#endif
	} else { // remove
	  fset.erase(-frag);
#ifdef DEBUG  
	  Rprintf("removed frag %d, size=%d\n",-frag,fset.size());
#endif
	}
	j++;
      }
#ifdef DEBUG  
      Rprintf("i=%d j=%d\n",i,j);
#endif
      if(return_list) {
	if(fset.empty() || (return_unique && fset.size()>1)) {
	  // assign null list?
	} else {
	  SEXP fil_R;
	  PROTECT(fil_R=allocVector(INTSXP,fset.size()));  np++;
	  int* fil=INTEGER(fil_R);
	  int k=0;
	  for(set<int>::const_iterator ki=fset.begin();ki!=fset.end();++ki) {
	    fil[k]=*ki; k++;
	  }
	  SET_VECTOR_ELT(nv, i, fil_R);
	  UNPROTECT(1); np--;
	}
      } else {
	if(return_point_counts) {
	  for(set<int>::const_iterator ki=fset.begin();ki!=fset.end();++ki) {
	    i_nv[*ki-1]++;
	  }
	} else {
	  if(fset.empty() || (return_unique && fset.size()>1)) {
	    i_nv[i]=-1;
	  } else {
	    i_nv[i]=*fset.begin();
	  }
	}
      }
    }

    UNPROTECT(np);
    return(nv);
  }


  SEXP expuni_lr(SEXP x_R,      // positions and their number (assumed sorted in ascending order)
		 SEXP mdist_R,  // max distance at which points should be considered
		 SEXP lambda_R,  // lambda value
		 SEXP spos_R,  // starting position
		 SEXP epos_R,  // ending position
		 SEXP step_R,  // step size
		 SEXP return_peaks_R, // whether peak positions should be returned, or entire score vector
		 SEXP min_peak_lr_R // min peak height (lr)
		 ) 
  {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    int* x=INTEGER(x_R);
    int nx=LENGTH(x_R);
    int mdist=INTEGER(mdist_R)[0];
    double lambda=*(REAL(lambda_R));

    int return_peaks=*(INTEGER(return_peaks_R));
    double min_peak=*(REAL(min_peak_lr_R));

    int spos=*(INTEGER(spos_R));
    int epos=*(INTEGER(epos_R));
    int step=*(INTEGER(step_R));

    int nsteps=(int) (epos-spos)/step;


#ifdef DEBUG  
    Rprintf("n=%d; lambda=%f; mdist=%d; spos=%d; epos=%d; step=%d; nsteps=%d\n",nx,lambda,mdist,spos,epos,step,nsteps);
#endif

    
    SEXP nv;
    double *d_nv;
    if(!return_peaks) {
      PROTECT(nv=allocVector(REALSXP,nsteps+1)); 
      d_nv=REAL(nv);
    }


    int i=0; // current index of the first point being used in the calculations
    int j=0; // current index of the last point being used in the calculations
    int sx=0; // current sum of all positions
    int n=0;
    
    for(int k=0; k<=nsteps; k++) {
      int cpos=spos+k*step;
      // increase i until x[i]>=cpos-mdist; remove x from sx; decrement n;
      while(i<nx && x[i]<(cpos-mdist)) {
	n--; sx-=x[i]; i++;
	//Rprintf("incremented i: i=%d; n=%d; sx=%d; cpos-mdist=%d; x[i-1]=%d\n",i,n,sx,cpos-mdist,x[i-1]);
      }
      //Rprintf("stable i: i=%d; n=%d; sx=%d; cpos-mdist=%d; x[i-1]=%d\n",i,n,sx,cpos-mdist,x[i-1]);

      //if(i>j) { j=i; }

      // increase j until x[j]>cpos
      while(j<nx && x[j]<=cpos) {
	n++; sx+=x[j]; j++;
	//Rprintf("incremented j: j=%d; n=%d; sx=%d; cpos=%d; x[j-1]=%d\n",j,n,sx,cpos,x[j-1]);
      }
      //Rprintf("stable j: j=%d; n=%d; sx=%d; cpos=%d; x[j-1]=%d\n",j,n,sx,cpos,x[j]);
      
      // calculate lr
      d_nv[k]=((double)(1-n))*log(lambda)-lambda*((double)(n*(cpos+1)-sx));
      //Rprintf("recorded lr[%d]=%f\n",k-1,d_nv[k-1]);
    }
    UNPROTECT(1);
    return(nv);
  }


  SEXP allpdist(SEXP x_R,SEXP max_dist_R) {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    double* x=REAL(x_R);
    int nx=LENGTH(x_R);
    double max_dist=*REAL(max_dist_R);
#ifdef DEBUG  
    Rprintf("n=%d; max_dist=%d\n",nx,max_dist);
#endif

    vector<double> dist;
    
    for(int i=0;i<nx;i++) {
      for(int j=i+1;j<nx;j++) {

	double d=x[j]-x[i];
#ifdef DEBUG  
	Rprintf("i=%d; j=%d; d=%f\n",i,j,d);
#endif
	if(d<=max_dist) {
	  dist.push_back(d);
	} else {
	  break;
	}
      }
    }

    SEXP nv;
    PROTECT(nv=allocVector(REALSXP,dist.size())); 
    double* i_nv=REAL(nv);
    int i=0;
    for(vector<double> ::const_iterator pi=dist.begin();pi!=dist.end();++pi) {
      i_nv[i++]=*pi;
    }
  
    UNPROTECT(1);
    return(nv);
  }

  // same as above, but for two different sets
  SEXP allxpdist(SEXP x_R,SEXP y_R, SEXP max_dist_R) {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    double* x=REAL(x_R);
    double* y=REAL(y_R);
    int nx=LENGTH(x_R);
    int ny=LENGTH(y_R);
    double max_dist=*REAL(max_dist_R);
#ifdef DEBUG  
    Rprintf("nx=%d; ny=%d; max_dist=%d\n",nx,ny,max_dist);
#endif

    vector<double> dist;
    int yi=0; // latest y start index

    for(int i=0;i<nx;i++) {
      // adjust yi so that yi>=x[i]-max_dist_R
      while(y[yi]<(x[i]-max_dist) && yi<ny) { yi++; }
      if(yi==ny) { break; }

      for(int j=yi;j<ny;j++) {
        double d=y[j]-x[i];
#ifdef DEBUG  
        Rprintf("i=%d; j=%d; d=%f\n",i,j,d);
#endif
        if(d<=max_dist) {
          dist.push_back(d);
        } else {
          break;
        }
      }
    }

    SEXP nv;
    PROTECT(nv=allocVector(REALSXP,dist.size()));
    double* i_nv=REAL(nv);
    int i=0;
    for(vector<double> ::const_iterator pi=dist.begin();pi!=dist.end();++pi) {
      i_nv[i++]=*pi;
    }

    UNPROTECT(1);
    return(nv);
  }

  // returns a vector giving for each point,
  // number of points within a given max_dist
  SEXP nwithindist(SEXP x_R,SEXP max_dist_R) {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    double* x=REAL(x_R);
    int nx=LENGTH(x_R);
    double max_dist=*REAL(max_dist_R);

    SEXP nv;
    PROTECT(nv=allocVector(REALSXP,nx)); 
    double* i_nv=REAL(nv);
    for(int i=0;i<nx;i++) { i_nv[i]=0; }

#ifdef DEBUG  
    Rprintf("n=%d; max_dist=%d\n",nx,max_dist);
#endif

    for(int i=0;i<nx;i++) {
      for(int j=i+1;j<nx;j++) {

	double d=x[j]-x[i];
#ifdef DEBUG  
	Rprintf("i=%d; j=%d; d=%f\n",i,j,d);
#endif
	if(d<=max_dist) {
	  i_nv[i]++;
	  i_nv[j]++;
	} else {
	  break;
	}
      }
    }
  
    UNPROTECT(1);
    return(nv);
  }




  // given a list of sorted chromosome signal and background vectors (unscaled), determine 
  // cluster contigs exceeding thr poisson P value, based on a whs window size,
  // and satisfying mcs cluster size
  SEXP find_poisson_enrichment_clusters(SEXP pos_R,SEXP flag_R,SEXP wsize_R,SEXP thr_R,SEXP mcs_R,SEXP bgm_R,SEXP mintag_R,SEXP either_R) {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    double* pos=REAL(pos_R);
    int* flag=INTEGER(flag_R);
    int nt=LENGTH(pos_R);
    
    int mcs=*INTEGER(mcs_R);
    int wsize=*INTEGER(wsize_R);
    int either=*INTEGER(either_R);
    double thr=REAL(thr_R)[0];
    double bgm=REAL(bgm_R)[0];
    double mintag=REAL(mintag_R)[0];

#ifdef DEBUG  
    Rprintf("nt=%d; wsize=%d; thr=%f; mcs=%d; min.tag=%f; bgm=%f\n",nt,wsize,thr,mcs,mintag,bgm);
#endif
    
    
    vector< pair<double,double> > contigs;
    
    // running indecies (start and end)
    int si=0; 
    int ei=0;
    
    // current window coordinate
    double ws=pos[0];
    
    // current window tag counts
    int cc[2]={0,0};


    if(nt>0) {
      cc[flag[si]]++;    
      // increment window end
      while(ei<(nt-1) && (pos[ei+1]-ws) <= wsize) {
	ei++;
	cc[flag[ei]]++;
      }


      // cluster start,end positions
      double cs,ce;
      int inclust=0;

      while(si<nt-1) {
      
	if((pos[si+1]-ws) > (pos[ei+1] - ws - wsize) && ei!=(nt-1)) {
	  // move end boudnary
	  ei++;
	  ws=pos[ei]-wsize;
	  cc[flag[ei]]++;
	  while(ei<(nt-1) && pos[ei+1]==ws+wsize) {
	    ei++;
	    cc[flag[ei]]++;
	  }
	
	  // increment window start
	  while(si<(nt-1) && pos[si] < ws) {
	    cc[flag[si]]--;
	    si++;
	  }

	} else {
	  // move up start boundary
	  ws=pos[si+1];
	  cc[flag[si]]--;
	  si++;
	  while(si<(nt-1) && pos[si+1]==ws) {
	    cc[flag[si]]--;
	    si++;
	  }
	
	  // increment window end
	  while(ei<(nt-1) && (pos[ei+1] - ws) <= wsize) {
	    ei++;
	    cc[flag[ei]]++;
	  }

	}

	// calculate z score
	double dc0=((double)cc[0])+0.5;
	double dc1=((double)cc[1])+0.5;
	double rte=dc0+dc1-0.25*thr*thr;
	double lb;
	if(rte<=0) { 
	  lb=0; 
	} else {
	  lb=(sqrt(dc1*dc0) - 0.5*thr*sqrt(rte))/(dc0 - 0.25*thr*thr);
	  if(lb<0) { lb=0; }
	  lb*=lb;
	}

	//Rprintf("%f=f(%f,%f,%f); %f=f(%f,%f,%f)\n",lb,1.0-thr,2.0*dc1,2.0*dc0,ub,thr,2.0*dc1,2.0*dc0);
      
#ifdef DEBUG  
	//double ub=gsl_cdf_fdist_Qinv(thr,2.0*dc1,2.0*dc0)*dc1/dc0;
	double ub=(sqrt(dc1*dc0) + 0.5*thr*sqrt(rte))/(dc0 - 0.25*thr*thr);
	ub*=ub;
	Rprintf("s=%d (%f); e=%d (%f); window: %f-%f; cc=[%d,%d]; lb=%f; ub=%f\n",si,pos[si],ei,pos[ei],ws,ws+wsize,cc[0],cc[1],lb,ub);
#endif
      
	int bc=lb>=bgm && cc[1]>=mintag;
	if(either) {
	  bc=lb>=bgm || cc[1]>=mintag;
	}
	if(bc) {
	  if(inclust) {
	    double nce=ws+wsize/2.0;
	    if(nce-ce > wsize/2.0) {
	      // next point is too far removed, end cluster
	      if(ce-cs >= mcs) {
		contigs.push_back(pair<double,double>(cs,ce));
#ifdef DEBUG  
		Rprintf("recorded cluster %f-%f\n",cs,ce);
#endif
	      }
	      inclust=0;
	    } else {
	      ce=nce;
	    }
	  } else {
	    inclust=1;
	    cs=ws+wsize/2.0;
	    ce=cs;
	  }
	} else {
	  if(inclust) {
	    if(ce-cs >= mcs) {
	      contigs.push_back(pair<double,double>(cs,ce));
#ifdef DEBUG  
	      Rprintf("recorded cluster %f-%f\n",cs,ce);
#endif
	    }
	    inclust=0;
	  }
	}
    
      }

      if(inclust) {
	if(ce-cs >= mcs) {
	  contigs.push_back(pair<double,double>(cs,ce));
#ifdef DEBUG  
	  Rprintf("recorded cluster %f-%f\n",cs,ce);
#endif
	}
	inclust=0;
      }
    }
    
    SEXP cs_R,ce_R;
    PROTECT(cs_R=allocVector(REALSXP,contigs.size())); 
    PROTECT(ce_R=allocVector(REALSXP,contigs.size())); 
    double* csa=REAL(cs_R);
    double* cea=REAL(ce_R);

    int i=0;
    for(vector< pair<double,double> >::const_iterator ci=contigs.begin(); ci!=contigs.end();++ci) {
      csa[i]=ci->first;
      cea[i]=ci->second;
      i++;
    }
    
    SEXP ans_R, names_R;
    PROTECT(names_R = allocVector(STRSXP, 2));
    SET_STRING_ELT(names_R, 0, mkChar("s"));
    SET_STRING_ELT(names_R, 1, mkChar("e"));
    
    PROTECT(ans_R = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans_R, 0, cs_R);
    SET_VECTOR_ELT(ans_R, 1, ce_R);
    setAttrib(ans_R, R_NamesSymbol, names_R);
  
    UNPROTECT(4);
    return(ans_R);

  }


  // finds intersection between a list of regions
  // the flag has +n/-n value, corresponding to the start/end of a segment in n-th regionset
  // max_val: 1 - report max overlapping value, -1: report min, 0 - don't look at values
  // returns: $s, $e, ($v) lists
  SEXP region_intersection(SEXP n_R,SEXP pos_R,SEXP flags_R,SEXP vals_R,SEXP max_val_R,SEXP union_R) {
    const int max_val=*INTEGER(max_val_R);
    const int unionr=*INTEGER(union_R);
    const int n=*INTEGER(n_R);
    double* pos=REAL(pos_R);
    int* flags=INTEGER(flags_R);
    double* val=REAL(vals_R);
    
#ifdef DEBUG  
    Rprintf("n=%d; npos=%d; max_val=%d\n",n,LENGTH(pos_R),max_val);
#endif

    int s[n]; // flag status for each set
    double mv[n]; // max/min value of current clusters

    for(int i=0;i<n;i++) { s[i]=0; }
    
    vector<double> starts;
    vector<double> ends;
    vector<double> values;

    int start=-1;
    double mval=0;
    for(int i=0;i<LENGTH(pos_R);i++) {
      // update flags
      int f=flags[i];
      if(f>0) {
	s[abs(f)-1]++;
      } else {
	s[abs(f)-1]--;
      }
      
      if(max_val!=0 && val[i]*max_val > mval*max_val) { mval=val[i]; }

      // joined status
      int all;
      if(unionr) {
	all=0;
	for(int j=0;j<n;j++) { if(s[j]>0) { all=1; break;} }
      } else {
	all=1;
	for(int j=0;j<n;j++) { all=all & (s[j]>0); }
      }
      
      
      //Rprintf("i=%d; s=[",i);
      //for(int j=0;j<n;j++) { Rprintf("%d",s[j]); }
      //Rprintf("]; all=%d; start=%d\n",all,start);

      if(start>=0) {
	// in fragment
	if(!all) { 
	  // end fragment
	  starts.push_back(pos[start]);
	  ends.push_back(pos[i]);
	  start=-1;
	  if(max_val!=0) { values.push_back(mval); }

#ifdef DEBUG  
	  Rprintf("recorded new fragment (s=%f,e=%f,v=%f);\n",pos[start],pos[i],mval);
#endif
	}
      } else {
	// should a fragment be started?
	if(all) {
	  start=i;
	  if(max_val!=0) { mval=val[i]; }
#ifdef DEBUG  
	  Rprintf("starting new fragment (s=%f,i=%d);\n",pos[start],i);
#endif
	}
      }
    }
    SEXP cs_R,ce_R,cv_R;
    PROTECT(cs_R=allocVector(REALSXP,starts.size())); 
    PROTECT(ce_R=allocVector(REALSXP,ends.size())); 
    
    double* csa=REAL(cs_R);
    int i=0;
    for(vector<double>::const_iterator ci=starts.begin(); ci!=starts.end(); ++ci) {
      csa[i]=*ci; i++;
    }

    csa=REAL(ce_R);
    i=0;
    for(vector<double>::const_iterator ci=ends.begin(); ci!=ends.end(); ++ci) {
      csa[i]=*ci; i++;
    }
    
    if(max_val!=0) {
      PROTECT(cv_R=allocVector(REALSXP,values.size())); 
      csa=REAL(cv_R);
      i=0;
      for(vector<double>::const_iterator ci=values.begin(); ci!=values.end(); ++ci) {
	csa[i]=*ci; i++;
      }
    }

    SEXP ans_R, names_R;
    if(max_val!=0) {
      PROTECT(names_R = allocVector(STRSXP, 3));
      SET_STRING_ELT(names_R, 0, mkChar("s"));
      SET_STRING_ELT(names_R, 1, mkChar("e"));
      SET_STRING_ELT(names_R, 2, mkChar("v"));

      PROTECT(ans_R = allocVector(VECSXP, 3));
      SET_VECTOR_ELT(ans_R, 0, cs_R);
      SET_VECTOR_ELT(ans_R, 1, ce_R);
      SET_VECTOR_ELT(ans_R, 2, cv_R);
    } else {
      PROTECT(names_R = allocVector(STRSXP, 2));
      SET_STRING_ELT(names_R, 0, mkChar("s"));
      SET_STRING_ELT(names_R, 1, mkChar("e"));

      PROTECT(ans_R = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(ans_R, 0, cs_R);
      SET_VECTOR_ELT(ans_R, 1, ce_R);
    }
    
    setAttrib(ans_R, R_NamesSymbol, names_R);
    
    if(max_val!=0) {
      UNPROTECT(5);
    } else {
      UNPROTECT(4);
    }
    return(ans_R);
  }

}

