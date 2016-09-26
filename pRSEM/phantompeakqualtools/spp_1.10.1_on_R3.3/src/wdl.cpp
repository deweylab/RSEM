#include <vector>
#include <string.h>
#include <iostream>
#include <string>
#include <set>

extern "C" {
// pliu 20160911
//#include "R.h"
//#include "Rmath.h"
#include <math.h>
//////
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 

//#define DEBUG 1

extern "C" {

  /************************************************************************/
  /*
   * lwcc - calculate local window cross-correlation
   */

  SEXP lwcc(SEXP x_R, // positive strand hist 
	    SEXP y_R, // negative strand hist of the same length
	    SEXP osize_R,       // outer boundary distance
	    SEXP isize_R,        // inner boundary distance
	    SEXP return_peaks_R, // whether all correlation values, or just peaks should be returned
	    SEXP min_peak_dist_R, // distance between closest peaks
	    SEXP min_peak_val_R, // min peak threshold
	    SEXP tag_weight_R,  // tag weight
	    SEXP bg_subtract_R, // a flag whether do background subtractio
	    SEXP bgp_R, // optional background hist for positive strand
	    SEXP bgn_R, // optional background hist for negative strand
	    SEXP bg_wsize_R, // window size for the background counts
	    SEXP bg_weight_R, // optional weighting for the background tags, must compensate for window size difference (including is cutout)
	    SEXP round_up_R // whether to round up fractional signal tag counts
	    )
  {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    
    int is=INTEGER(isize_R)[0];
    int os=INTEGER(osize_R)[0];
    double rs=((double)(2*os+1));
    int* x=INTEGER(x_R);
    int* y=INTEGER(y_R);
    int n_x=LENGTH(x_R);

    // background-related 
    int* bgp=INTEGER(bgp_R);
    int* bgn=INTEGER(bgn_R);
    int bg_whs=INTEGER(bg_wsize_R)[0];

    int return_peaks=*(INTEGER(return_peaks_R));
    double min_peak_val=*(REAL(min_peak_val_R));
    int min_peak_dist=*(INTEGER(min_peak_dist_R));
    double tag_weight=*(REAL(tag_weight_R));

    const int round_up=*(INTEGER(round_up_R));
    const int bg_subtract=*(INTEGER(bg_subtract_R));
    const double bg_weight=*(REAL(bg_weight_R));

    int i; // point at which the value is being calculated
    int start=os;
    int end=n_x-os-1;

    // bg tag counts within bg window
    int bg_pn1=0;
    int bg_nn1=0;
    int bg_pn2=0;
    int bg_nn2=0;


  
    // illustration for counting:
    //
    // 012345678901234567890123456789012
    // ==========------|------==========
    //
    //  osize=16; isize=6; 


    SEXP nv;
    double *d_nv;
    vector<int> ppos;
    vector<double> pval;
    if(!return_peaks) {
      PROTECT(nv=allocVector(REALSXP,n_x)); 
      d_nv=REAL(nv);
      for(int i=0;i<n_x;i++) {
	d_nv[i]=0;
      }
    }

#ifdef DEBUG  
    Rprintf("start=%d end=%d tag_weight=%f\n", start,end,tag_weight);
    Rprintf("x[1]=%d x[2]=%d y[1]=%d y[2]=%d\n",x[1],x[2],y[1],y[2]);
#endif

    int lpp=-1; // last peak position
    double lpv=-1e3; // last peak value
    
    double ppv=-1e3; // last value
    double pppv=-11e-3; // value before last

    int pn1,pn2,nn1,nn2;

    
    if(bg_subtract) {
      // pre-initialize background tag counts, 
      for(int i=0;i<bg_whs;i++) {
	if(i<n_x) {
	  bg_pn2+=bgp[i];
	  bg_nn2+=bgn[i];
	}
      }
    }


    for(i=0;i<end;i++) {
#ifdef DEBUG  
      //Rprintf("i=%d ", i);
#endif
      
      if(bg_subtract) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

      if(i >= start) {
	// update counts, taking into account masked out regions
	pn1=pn2=nn1=nn2=0;
	
	for(int k=0;k<=(os-is);k++) {
	  int xp1=x[i-os+k];
	  int xp2=x[i+os-k];
	  int xn1=y[i+os-k];
	  int xn2=y[i-os+k];

	  if(xp1!=-1 && xn1!=-1) {
	    pn1+=xp1;
	    nn1+=xn1;
	  }
	  if(xp2!=-1 && xn2!=-1) {
	    pn2+=xp2;
	    nn2+=xn2;
	  }
	}
      
	// calculate the means
	double mp=((double)(pn1+pn2))/rs;
	double mn=((double)(pn1+pn2))/rs;
#ifdef DEBUG  
	Rprintf("mp=%f mn=%f\n",mp,mn);
#endif
	// calculate correlation
	double varp=0;
	double varn=0;
	double num=0;
	double val=-1e3;
	if(mp>0 & mn>0) {
	  for(int k=0;k<=(os-is);k++) {
	    int xp1=x[i-os+k];
	    int xp2=x[i+os-k];
	    int xn1=y[i+os-k];
	    int xn2=y[i-os+k];

	    
	    if(xp1!=-1 && xn1!=-1) {  
	      double nnp1=((double) xp1)-mp;
	      double nnn1=((double) xn1)-mn;
	      num+=nnp1*nnn1;
	      varp+=nnp1*nnp1;
	      varn+=nnn1*nnn1;
	    }
	    
	    if(xp2!=-1 && xn2!=-1) {
	      double nnp2=((double) xp2)-mp;
	      double nnn2=((double) xn2)-mn;
	      num+=nnp2*nnn2;
	      varp+=nnp2*nnp2;
	      varn+=nnn2*nnn2;
	    }

	  }
	  double tagw;
	  double spn1=((double)pn1)*tag_weight;
	  double snn1=((double)nn1)*tag_weight;
	  double spn2=((double)pn2)*tag_weight;
	  double snn2=((double)nn2)*tag_weight;
	  if(round_up) {
	    if(pn1>0 && spn1<1) { spn1=1.0; }
	    //if(pn2>0 && spn2<1) { spn2=1.0; }
	    if(nn1>0 && snn1<1) { snn1=1.0; }
	    //if(nn2>0 && snn2<1) { snn2=1.0; }
	  }

	  if(bg_subtract) {
	    spn1-=((double)bg_pn1)*bg_weight;
	    snn1-=((double)bg_nn2)*bg_weight;
	    spn2-=((double)bg_pn2)*bg_weight;
	    snn2-=((double)bg_nn1)*bg_weight;

	    if(spn2<0) spn2=0;
	    if(snn2<0) snn2=0;
	    
	    if(spn1>0 && snn1>0) {
	      tagw=(2.0*sqrt(spn1*snn1)-(spn2+snn2+1.0));
	    } else {
	      tagw=-(spn2+snn2+1.0);
	    }
	    //cout<<"bg_pn1="<<bg_pn1<<"; bg_pn2="<<bg_pn2<<"; bg_nn1="<<bg_nn1<<"; bg_nn2="<<bg_nn2<<endl;
	  } else {
	    tagw=2.0*sqrt(spn1*snn1)-(spn2+snn2);
	  }

	  if(tagw<0) {
	    val=0.0; 
	  } else {
	    if(num==0.0) {
	      val=0;
	    } else {
	      val=num/(sqrt(varp*varn));
	    }
	    val=val*sqrt(tagw) + tagw;

	  }
	  //cout<<"val="<<val<<endl;

#ifdef DEBUG  
        Rprintf("pn1=%d pn2=%d nn1=%d nn2=%d tag.weight=%f tagw=%f\n",pn1,pn2,nn1,nn2,tag_weight,tagw);
	Rprintf("tagw=%f varp=%f varn=%f num=%f cor=%f val=%f\n",tagw,varp,varn,num,num/sqrt(varp*varn),val);
#endif
	}


	
	if(return_peaks) {
	  // determine if previous position was a peak
	  if(ppv>min_peak_val && ppv>val && ppv>pppv) {
	    if(lpp>0 && (i-lpp+1)>min_peak_dist) {
	      // record previous peak position
	      ppos.push_back(lpp);
	      pval.push_back(lpv);
#ifdef DEBUG  
	      Rprintf("recording peak x=%d y=%f d=%d\n",lpp,lpv,(i-lpp));
#endif	    
	      lpp=i-1; lpv=ppv;
#ifdef DEBUG  
	      Rprintf("updated peak to x=%d y=%f\n",lpp,lpv);
#endif	    
	    } else {
	      if(ppv>lpv) {
		// update last peak positions
#ifdef DEBUG  
		Rprintf("skipping peak x=%d y=%f d=%d in favor of x=%d y=%f\n",lpp,lpv,(i-lpp),i-1,ppv);
#endif
		lpp=i-1; lpv=ppv;
	      }
	    }
	  }

	  // update previous values
	  if(val!=ppv) {
	    pppv=ppv; ppv=val;
	  }
	} else {
	  d_nv[i]=val;
	}
      }
    }

    if(return_peaks) {
      // record last position
      if(lpp>0) {
#ifdef DEBUG  
	Rprintf("recording last peak x=%d y=%f\n",lpp,lpv);
#endif
	ppos.push_back(lpp);
	pval.push_back(lpv);
      }

      SEXP rpp_R,rpv_R;
      PROTECT(rpp_R=allocVector(INTSXP,ppos.size())); 
      PROTECT(rpv_R=allocVector(REALSXP,ppos.size())); 
      int* rpp=INTEGER(rpp_R);
      double* rpv=REAL(rpv_R);

      for(int i=0;i<ppos.size();i++) {
	rpp[i]=ppos[i];
	rpv[i]=pval[i];
      }
    
      SEXP ans_R, names_R;
      PROTECT(names_R = allocVector(STRSXP, 2));
      SET_STRING_ELT(names_R, 0, mkChar("x"));
      SET_STRING_ELT(names_R, 1, mkChar("v"));
    
      PROTECT(ans_R = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(ans_R, 0, rpp_R);
      SET_VECTOR_ELT(ans_R, 1, rpv_R);
      setAttrib(ans_R, R_NamesSymbol, names_R);
  
      UNPROTECT(4);
      return(ans_R);
    } else {
      UNPROTECT(1);
      return(nv);
    }

  }



  /************************************************************************/
  /*
   * wtd - window tag difference implementation
   */

  SEXP wtd(SEXP x_R, // positive strand hist 
	   SEXP y_R, // negative strand hist of the same length
	   SEXP wsize_R,       // outer boundary distance
	   SEXP return_peaks_R, // whether all correlation values, or just peaks should be returned
	   SEXP min_peak_dist_R, // distance between closest peaks
	   SEXP min_peak_val_R, // min peak threshold
	   SEXP direct_count_R, // whether tag weighting should not be done
	   SEXP tag_weight_R,  // tag weight
	   SEXP ignore_masking_R,  // whether to ignore masked regions
	   SEXP bg_subtract_R, // a flag whether do background subtractio
	   SEXP bgp_R, // optional background hist for positive strand
	   SEXP bgn_R, // optional background hist for negative strand
	   SEXP bg_wsize_R, // window size for the background counts
	   SEXP bg_weight_R, // optional weighting for the background tags, must compensate for window size difference
	   SEXP round_up_R // whether to round up fractional signal tag counts
	   )
  {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    
    int whs=INTEGER(wsize_R)[0];
    int* x=INTEGER(x_R);
    int* y=INTEGER(y_R);
    int n_x=LENGTH(x_R);
    
    // background-related 
    int* bgp=INTEGER(bgp_R);
    int* bgn=INTEGER(bgn_R);
    int bg_whs=INTEGER(bg_wsize_R)[0];
    

    const int return_peaks=*(INTEGER(return_peaks_R));
    const int direct_count=*(INTEGER(direct_count_R));
    const int ignore_masking=*(INTEGER(ignore_masking_R));
    const double min_peak_val=*(REAL(min_peak_val_R));
    const int min_peak_dist=*(INTEGER(min_peak_dist_R));
    const double tag_weight=*(REAL(tag_weight_R));
    
    const int round_up=*(INTEGER(round_up_R));
    const int bg_subtract=*(INTEGER(bg_subtract_R));
    const double bg_weight=*(REAL(bg_weight_R));
    
    int i; // point at which the value is being calculated
    int start=whs+1;
    int end=n_x-whs-1;

    // tag counts to calculate the means
    int pn1=0;
    int pn2=0;
    int nn1=0;
    int nn2=0;

    // bg tag counts within bg window
    int bg_pn1=0;
    int bg_pn2=0;
    int bg_nn1=0;
    int bg_nn2=0;
    
    SEXP nv;
    double *d_nv;
    vector<int> ppos;
    vector<double> pval;
    if(!return_peaks) {
      PROTECT(nv=allocVector(REALSXP,n_x)); 
      d_nv=REAL(nv);
      for(int i=0;i<n_x;i++) {
	d_nv[i]=0;
      }
    }

#ifdef DEBUG  
    Rprintf("whs=%d start=%d end=%d tag_weight=%f ignore_masing=%d\n", whs, start,end,tag_weight,ignore_masking);
    Rprintf("x[1]=%d x[2]=%d y[1]=%d y[2]=%d\n",x[1],x[2],y[1],y[2]);
#endif

    int lpp=-1; // last peak position
    double lpv=-1000; // last peak value
    
    double ppv=-1000; // last value
    int ppl=-1; // position of the last value
    double pppv=-1000; // value before last


    if(ignore_masking==1) {
      for(int i=0;i<whs;i++) {
	pn1+=x[i];
	pn2+=x[i+whs+1];
	nn1+=y[i];
	nn2+=y[i+whs+1];

      }
    }

    if(bg_subtract) {
      // pre-initialize background tag counts, 
      for(int i=0;i<bg_whs;i++) {
	if(i<n_x) {
	  bg_pn2+=bgp[i];
	  bg_nn2+=bgn[i];
	}
      }
      // increment center of background count window to the start position
      for(int i=0;i<start;i++) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

    }

    
#ifdef DEBUG  
    Rprintf("initialization: i=%d pn1=%d, pn2=%d, nn1=%d, nn2=%d", i,pn1,pn2,nn1,nn2);
#endif

    for(i=start;i<end;i++) {
      if(bg_subtract) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

      // update counts
      if(ignore_masking==1) {
	pn1+=x[i-1]-x[i-whs-1];
	pn2+=x[i+whs]-x[i-1];
	nn1+=y[i-1]-y[i-whs-1];
	nn2+=y[i+whs]-y[i-1];

      } else {

	pn1=pn2=nn1=nn2=0;
	
	for(int k=0;k<whs;k++) {
	  int xp1=x[i-k-1];
	  int xp2=x[i+k];
	  int xn1=y[i-k-1];
	  int xn2=y[i+k];

	  // omit masked positions
	  if(xp1!=-1 && xn1!=-1 && xp2!=-1 && xn2!=-1) {
	    pn1+=xp1;
	    nn1+=xn1;
	    pn2+=xp2;
	    nn2+=xn2;
	  }
	}
      }

      double val;
      double spn1=((double)pn1)*tag_weight;
      double snn1=((double)nn1)*tag_weight;
      double spn2=((double)pn2)*tag_weight;
      double snn2=((double)nn2)*tag_weight;
      if(round_up) {
	if(pn1>0 && spn1<1) { spn1=1.0; }
	//if(pn2>0 && spn2<1) { spn2=1.0; }
	//if(nn1>0 && snn1<1) { snn1=1.0; }
	if(nn2>0 && snn2<1) { snn2=1.0; }
      }

      if(direct_count) {
	val=spn1+snn2;
	if(round_up && val<1) {
	  val=1.0;
	}
	if(bg_subtract) {
	  val-=((double) (bg_pn1+bg_nn2))*bg_weight;
	}
      } else {
	if(bg_subtract) {
	  spn1-=((double)bg_pn1)*bg_weight;
	  snn1-=((double)bg_nn1)*bg_weight;
	  spn2-=((double)bg_pn2)*bg_weight;
	  snn2-=((double)bg_nn2)*bg_weight;

	  if(spn2<0) spn2=0;
	  if(snn1<0) snn1=0;

	  if(spn1>0 && snn2>0) {
	    val=(2.0*sqrt(spn1*snn2)-(spn2+snn1+1.0));
	  } else {
	    val=-(spn2+snn1+1.0);
	  }
	} else {
	  val=2.0*sqrt(spn1*snn2)-(spn2+snn1+tag_weight);
	}
      }	
      //double val=sqrt(pn1*nn2);
      //if(pn2>nn1) { val-=pn2; } else { val-=pn1; }
#ifdef DEBUG  
      Rprintf("update: i=%d pn1=%d pn2=%d nn1=%d nn2=%d val=%f\n",i,pn1,pn2,nn1,nn2,val);
#endif
      
      if(return_peaks) {
	// determine if previous position was a peak
	if(ppv>min_peak_val && ppv>val && ppv>pppv) {
	  if(lpp>0 && (i-lpp+1)>min_peak_dist) {
	    // record previous peak position
	    ppos.push_back(lpp);
	    pval.push_back(lpv);
#ifdef DEBUG  
	    Rprintf("recording peak x=%d y=%f d=%d\n",lpp,lpv,(i-lpp));
#endif	    
	    if(ppl!=-1 && ppl!=i-1) {
	      lpp=(int) round((ppl+i-1)/2);
	    } else {
	      lpp=i-1;
	    }
	    lpv=ppv;
#ifdef DEBUG  
	    Rprintf("updated peak to x=%d y=%f\n",lpp,lpv);
#endif	    
	  } else {
	    if(ppv>lpv) {
	      // update last peak positions
#ifdef DEBUG  
	      Rprintf("skipping peak x=%d y=%f d=%d in favor of x=%d y=%f\n",lpp,lpv,(i-lpp),i-1,ppv);
#endif
	      if(ppl!=-1 && ppl!=i-1) {
		lpp=(int) round((ppl+i-1)/2);
	      } else {
		lpp=i-1;
	      }
	      lpv=ppv;
	    }
	  }
	}
	
	// update previous values
	if(val!=ppv) {
	  pppv=ppv; ppv=val; ppl=i;
	}
      } else {
	d_nv[i]=val;
      }
    }

    if(return_peaks) {
      // record last position
      if(lpp>0) {
#ifdef DEBUG  
	Rprintf("recording last peak x=%d y=%f\n",lpp,lpv);
#endif
	ppos.push_back(lpp);
	pval.push_back(lpv);
      }

      SEXP rpp_R,rpv_R;
      PROTECT(rpp_R=allocVector(INTSXP,ppos.size())); 
      PROTECT(rpv_R=allocVector(REALSXP,ppos.size())); 
      int* rpp=INTEGER(rpp_R);
      double* rpv=REAL(rpv_R);

      for(int i=0;i<ppos.size();i++) {
	rpp[i]=ppos[i];
	rpv[i]=pval[i];
      }
    
      SEXP ans_R, names_R;
      PROTECT(names_R = allocVector(STRSXP, 2));
      SET_STRING_ELT(names_R, 0, mkChar("x"));
      SET_STRING_ELT(names_R, 1, mkChar("v"));
    
      PROTECT(ans_R = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(ans_R, 0, rpp_R);
      SET_VECTOR_ELT(ans_R, 1, rpv_R);
      setAttrib(ans_R, R_NamesSymbol, names_R);
  
      UNPROTECT(4);
      return(ans_R);
    } else {
      UNPROTECT(1);
      return(nv);
    }

  }


}


