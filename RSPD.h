#ifndef RSPD_H_
#define RSPD_H_

#include<cstdio>
#include<cstring>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "Refs.h"
#include "simul.h"

const int RSPD_DEFAULT_B = 20;

class RSPD {
public:
	RSPD(bool estRSPD, int B = RSPD_DEFAULT_B) {
		this->estRSPD = estRSPD;
		this->B = B;

		pdf = new double[B + 2];
		cdf = new double[B + 2];

		//set initial parameters
		memset(pdf, 0, sizeof(double) * (B + 2)); // use B + 2 for evalCDF
		memset(cdf, 0, sizeof(double) * (B + 2));
		for (int i = 1; i <= B; i++) {
			pdf[i] = 1.0 / B;
			cdf[i] = i * 1.0 / B;
		}
	}

	~RSPD() {
		delete[] pdf;
		delete[] cdf;
	}

	RSPD& operator=(const RSPD& rv);

	void init();

	//fpos starts from 0
	void update(int fpos, int fullLen, double frac) {
		assert(estRSPD);

		if (fpos >= fullLen) return;  // if out of range, do not use this hit

		int i;
		double a = fpos * 1.0 / fullLen;
		double b;

		for (i = ((long long)fpos) * B / fullLen + 1; i < (((long long)fpos + 1) * B - 1) / fullLen + 1; i++) {
			b = i * 1.0 / B;
			pdf[i] += (b - a) * fullLen * frac;
			a = b;
		}
		b = (fpos + 1.0) / fullLen;
		pdf[i] += (b - a) * fullLen * frac;
	}

	void finish();

	double evalCDF(int fpos, int fullLen) {
		int i = ((long long)fpos) * B / fullLen;
		double val = fpos * 1.0 / fullLen * B;

		return cdf[i] + (val - i) * pdf[i + 1];
	}

	double getAdjustedProb(int fpos, int effL, int fullLen) {
		assert(fpos >= 0 && fpos < fullLen && effL <= fullLen);
		if (!estRSPD) return 1.0 / effL;
		double denom = evalCDF(effL, fullLen);
		return (denom >= EPSILON ? (evalCDF(fpos + 1, fullLen) - evalCDF(fpos, fullLen)) / denom : 0.0) ;
	}

	void collect(const RSPD&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation(int, Refs*);
	int simulate(simul*, int, int);
	void finishSimulation();

private:
	bool estRSPD;
	int B; // number of bins
	double *pdf, *cdf;

	int M;
	double **rspdDists;
};

RSPD& RSPD::operator=(const RSPD& rv) {
	if (this == &rv) return *this;
	if (B != rv.B) {
		delete[] pdf;
		delete[] cdf;
		pdf = new double[rv.B + 2];
		cdf = new double[rv.B + 2];
	}
	B = rv.B;
	memcpy(pdf, rv.pdf, sizeof(double) * (B + 2));
	memcpy(cdf, rv.cdf, sizeof(double) * (B + 2));

	return *this;
}

void RSPD::init() {
	assert(estRSPD);
	memset(pdf, 0, sizeof(double) * (B + 2));
	memset(cdf, 0, sizeof(double) * (B + 2));
}

void RSPD::finish() {
	double sum = 0.0;

	assert(estRSPD);

	for (int i = 1; i <= B; i++) {
		sum += pdf[i];
	}

	for (int i = 1; i <= B; i++) {
		pdf[i] /= sum;
		cdf[i] = cdf[i - 1] + pdf[i];
	}
}

void RSPD::collect(const RSPD& o) {
	assert(estRSPD);
	for (int i = 1; i <= B; i++) {
		pdf[i] += o.pdf[i];
	}
}

void RSPD::read(FILE *fi) {
	//release default space first
	delete[] pdf;
	delete[] cdf;

	int val;
	assert(fscanf(fi, "%d", &val) == 1);
	estRSPD = (val != 0);

	if (estRSPD) {
		assert(fscanf(fi, "%d", &B) == 1);
		pdf = new double[B + 2];
		cdf = new double[B + 2];
		memset(pdf, 0, sizeof(double) * (B + 2));
		memset(cdf, 0, sizeof(double) * (B + 2));
		for (int i = 1; i <= B; i++) {
			assert(fscanf(fi, "%lf", &pdf[i]) == 1);
			cdf[i] = cdf[i - 1] + pdf[i];
		}
	}
	else {
		B = RSPD_DEFAULT_B;
		pdf = new double[B + 2];
		cdf = new double[B + 2];
		memset(pdf, 0, sizeof(double) * (B + 2));
		memset(cdf, 0, sizeof(double) * (B + 2));
		for (int i = 1; i <= B; i++) {
			pdf[i] = 1.0 / B;
			cdf[i] = i * 1.0 / B;
		}
	}
}

void RSPD::write(FILE *fo) {
	fprintf(fo, "%d\n", estRSPD);
	if (estRSPD) {
		fprintf(fo, "%d\n", B);
		for (int i = 1; i < B; i++) {
			fprintf(fo, "%.10g ", pdf[i]);
		}
		fprintf(fo, "%.10g\n", pdf[B]);
	}
}

void RSPD::startSimulation(int M, Refs* refs) {
	if (!estRSPD) return;
	this->M = M;
	rspdDists = new double*[M + 1];
	rspdDists[0] = NULL;
	for (int i = 1; i <= M; i++) {
		int fullLen = refs->getRef(i).getFullLen();
		rspdDists[i] = new double[fullLen];
		memset(rspdDists[i], 0, sizeof(double) * fullLen);
		for (int j = 0; j < fullLen; j++) rspdDists[i][j] = evalCDF(j + 1, fullLen);
	}
}

int RSPD::simulate(simul *sampler, int sid, int effL) {
	if (estRSPD) return (rspdDists[sid][effL - 1] > 0.0 ? sampler->sample(rspdDists[sid], effL) : -1);
	return int(sampler->random() * effL);
}

void RSPD::finishSimulation() {
	if (!estRSPD) return;
	for (int i = 1; i <= M; i++) delete[] rspdDists[i];
	delete[] rspdDists;
}

#endif /* RSPD_H_ */
