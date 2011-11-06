#ifndef LENDIST_H_
#define LENDIST_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<algorithm>

#include "boost/math/distributions/normal.hpp"

#include "utils.h"
#include "simul.h"

class LenDist {
public:
	LenDist(int minL = 1, int maxL = 1000) {
		lb = minL - 1;
		ub = maxL;
		span = ub - lb;
		assert(span > 0);

		pdf = new double[span + 1];
		cdf = new double[span + 1];

		//set initial parameters
		pdf[0] = cdf[0] = 0.0;
		for (int i = 1; i <= span; i++) {
			pdf[i] = 1.0 / span;
			cdf[i] = i * 1.0 / span;
		}
	}

	~LenDist() {
		delete[] pdf;
		delete[] cdf;
	}

	LenDist& operator=(const LenDist&);

	void setAsNormal(double, double, int, int);

	void init();

	//the corresponding lb and ub are the original one
	void update(int len, double frac) {
		assert(len > lb && len <= ub);
		pdf[len - lb] += frac;
	}

	void finish();

	int getMinL() const { return lb + 1; }
	int getMaxL() const { return ub; }

	double getProb(int len) const {
		assert(len > lb && len <= ub);
		return pdf[len - lb];
	}

	//len : mate/fragment length
	//refL : reference sequence length, in fact, this is totLen for global length distribution
	double getAdjustedProb(int len, int refL) const {
		if (len <= lb || len > ub || refL <= lb) return 0.0;
		double denom = cdf[std::min(ub, refL) - lb];
		assert(denom >= EPSILON);
		return pdf[len - lb] / denom;
	}

	//len : length threshold, any length <= len should be calculated
	//refL : reference sequence length
	double getAdjustedCumulativeProb(int len, int refL) const {
		assert(len > lb && len <= ub && refL > lb);
		double denom = cdf[std::min(ub, refL) - lb];
		assert(denom >= EPSILON);
		return cdf[len - lb] / denom;
	}

	//for multi-thread usage
	void collect(const LenDist&);

	void read(FILE*);
	void write(FILE*);

	void copyTo(double*&, double*&, int&, int&, int&) const;
	
	int simulate(simul*, int);
	
 private:
	int lb, ub, span; // (lb, ub]
	double *pdf, *cdf;

	void trim();
};

LenDist& LenDist::operator=(const LenDist& rv) {
	if (this == &rv) return *this;
	if (span != rv.span) {
		delete[] pdf;
		delete[] cdf;
		pdf = new double[rv.span + 1];
		cdf = new double[rv.span + 1];
	}
	lb = rv.lb; ub = rv.ub; span = rv.span;
	memcpy(pdf, rv.pdf, sizeof(double) * (span + 1));
	memcpy(cdf, rv.cdf, sizeof(double) * (span + 1));

	return *this;
}

//Please give interger mean, thanks!
//minL: new minimum length, maxL: new maximum length
void LenDist::setAsNormal(double mean, double sd, int minL, int maxL) {
  int meanL = int(mean + .5); // assume meanL is a integer; if not, round to nearest number.    
  delete[] pdf;
  delete[] cdf;

  if (sd < EPSILON) {
    if (meanL < minL || meanL > maxL) {
      fprintf(stderr, "Length distribution's probability mass is not within the possible range! MeanL = %d, MinL = %d, MaxL = %d\n", meanL, minL, maxL);
      exit(-1);
    }
    span = 1;
    lb = meanL - 1; ub = meanL;
    pdf = new double[span + 1];
    cdf = new double[span + 1];
    pdf[0] = cdf[0] = 0.0;
    pdf[1] = cdf[1] = 1.0;

    return;
  }


  boost::math::normal norm(mean, sd);

  if (maxL - minL + 1 > RANGE) {
    if (meanL <= minL) maxL = minL + RANGE - 1;
    else if (meanL >= maxL) minL = maxL - RANGE + 1;
    else {
      double lg = mean - (minL - 0.5);
      double rg = (maxL + 0.5) - mean;
      double half = RANGE / 2.0;
      
      if (lg < half) { assert(rg > half); maxL = minL + RANGE - 1; }
      else if (rg < half) { assert(lg > half); minL = maxL - RANGE + 1; }
      else { minL = int(mean - half + 1.0); maxL = int(mean + half); }
    }
  }

  assert(maxL - minL + 1 <= RANGE);

  lb = minL - 1;
  ub = maxL;
  span = ub - lb;
  assert(span > 0);
  
  pdf = new double[span + 1];
  cdf = new double[span + 1];
  
  pdf[0] = cdf[0] = 0.0;
  
  double old_val, val, sum;
    
  sum = 0.0;
  old_val = boost::math::cdf(norm, minL - 0.5);
  for (int i = 1; i <= span; i++) {
    val = boost::math::cdf(norm, lb + i + 0.5);
    pdf[i] = val - old_val;
    sum += pdf[i];
    old_val = val;
  }
  assert(sum >= EPSILON);
  for (int i = 1; i <= span; i++) {
    pdf[i] /= sum;
    cdf[i] = cdf[i - 1] + pdf[i];
  }
  
  trim();
}

void LenDist::init() {
	memset(pdf, 0, sizeof(double) * (span + 1));
	memset(cdf, 0, sizeof(double) * (span + 1));
}

void LenDist::finish() {
	double sum = 0.0;

	for (int i = 1; i <= span; i++) {
		sum += pdf[i];
	}

	if (sum <= EPSILON) { fprintf(stderr, "No valid read to estimate the length distribution!\n"); exit(-1); }

	for (int i = 1; i <= span; i++) {
		pdf[i] = pdf[i] / sum;
		cdf[i] = cdf[i - 1] + pdf[i];
	}
	trim();
}


void LenDist::collect(const LenDist& o) {
	if (lb != o.lb || ub != o.ub) {
	  delete[] pdf;
	  delete[] cdf;
	  lb = o.lb; ub = o.ub; span = o.span;
	  pdf = new double[span + 1];
	  cdf = new double[span + 1];
	  memset(pdf, 0, sizeof(double) * (span + 1));
	  memset(cdf, 0, sizeof(double) * (span + 1));
	}
	for (int i = 1; i <= span; i++) {
		pdf[i] += o.pdf[i];
	}
}

void LenDist::read(FILE *fi) {
	//release default space first
	delete[] pdf;
	delete[] cdf;

	assert(fscanf(fi, "%d %d %d", &lb, &ub, &span) == 3);
	pdf = new double[span + 1];
	cdf = new double[span + 1];
	pdf[0] = cdf[0] = 0.0;
	for (int i = 1; i <= span; i++) {
	        assert(fscanf(fi, "%lf", &pdf[i]) == 1);
		cdf[i] = cdf[i - 1] + pdf[i];
	}

	trim();
}

void LenDist::write(FILE *fo) {
	fprintf(fo, "%d %d %d\n", lb, ub, span);
	for (int i = 1; i < span; i++) {
		fprintf(fo, "%.10g ", pdf[i]);
	}
	fprintf(fo, "%.10g\n", pdf[span]);
}

void LenDist::copyTo(double*& pdf, double*& cdf, int& lb, int& ub, int& span) const {
	lb = this->lb;
	ub = this->ub;
	span = this->span;

	pdf = new double[span + 1];
	memcpy(pdf, this->pdf, sizeof(double) * (span + 1));
	cdf = new double[span + 1];
	memcpy(cdf, this->cdf, sizeof(double) * (span + 1));
}

//refL = -1 means that this length is generated for noise isoform
int LenDist::simulate(simul* sampler, int refL) {
	int dlen;

	if (refL == -1) refL = ub;
	if (refL <= lb || cdf[(dlen = std::min(ub, refL) - lb)] <= 0.0) return -1;
	int len = lb + 1 + sampler->sample(cdf + 1, dlen);

	return len;
}

void LenDist::trim() {
  int newlb, newub;
  double *newpdf, *newcdf;

  for (newlb = 1; newlb <= span && pdf[newlb] < EPSILON; newlb++);
  newlb--;
  for (newub = span; newub > newlb && pdf[newub] < EPSILON; newub--);
  assert(newlb < newub);
  if (newlb == 0 && newub == span) return;

  span = newub - newlb;
  newpdf = new double[span + 1];
  memset(newpdf, 0, sizeof(double) * (span + 1));
  newcdf = new double[span + 1];
  memset(newcdf, 0, sizeof(double) * (span + 1));

  for (int i = 1; i <= span; i++) {
    newpdf[i] = pdf[i + newlb];
    newcdf[i] = cdf[i + newlb];
  }

  delete[] pdf;
  delete[] cdf;

  pdf = newpdf;
  cdf = newcdf;

  lb += newlb;
  ub = lb + span;
}

#endif /* LENDIST_H_ */
