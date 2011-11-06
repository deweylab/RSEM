#ifndef NOISEPROFILE_H_
#define NOISEPROFILE_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<string>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "simul.h"

class NoiseProfile {
public:
	NoiseProfile() {
		logp = 0.0;
		memset(c, 0, sizeof(c));
		memset(p, 0, sizeof(p));
	}

	NoiseProfile& operator=(const NoiseProfile&);

	void init();
	void updateC(const std::string&);
	void update(const std::string&, double frac);
	void finish();
	void calcInitParams();

	double getProb(const std::string&);
	double getLogP() { return logp; }

	void collect(const NoiseProfile&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int);
	void finishSimulation();

private:
	static const int NCODES = 5;

	double logp;
	double c[NCODES]; // counts in N0;
	double p[NCODES];

	double *pc; // for simulation
};

NoiseProfile& NoiseProfile::operator=(const NoiseProfile& rv) {
	if (this == &rv) return *this;
	logp = rv.logp;
	memcpy(c, rv.c, sizeof(rv.c));
	memcpy(p, rv.p, sizeof(rv.p));
	return *this;
}

void NoiseProfile::init() {
	memset(p, 0, sizeof(p));
}

void NoiseProfile::updateC(const std::string& readseq) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		++c[get_base_id(readseq[i])];
	}
}

void NoiseProfile::update(const std::string& readseq, double frac) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		p[get_base_id(readseq[i])] += frac;
	}
}

void NoiseProfile::finish() {
	double sum;

	logp = 0.0;
	sum = 0.0;
	for (int i = 0; i < NCODES; i++) sum += (p[i] + c[i]);
	if (sum <= EPSILON) return;
	for (int i = 0; i < NCODES; i++) {
		p[i] = (p[i] + c[i]) / sum;
		if (c[i] > 0.0) { logp += c[i] * log(p[i]); }
	}
}

void NoiseProfile::calcInitParams() {
	double sum;

	logp = 0.0;
	sum = 0.0;
	for (int i = 0; i < NCODES; i++) sum += (1.0 + c[i]);
	for (int i = 0; i < NCODES; i++) {
		p[i] = (1.0 + c[i]) / sum;
		if (c[i] > 0.0) { logp += c[i] * log(p[i]); }
	}
}

double NoiseProfile::getProb(const std::string& readseq) {
	double prob = 1.0;
	int len = readseq.size();

	for (int i = 0; i < len; i++) {
		prob *= p[get_base_id(readseq[i])];
	}

	return prob;
}

void NoiseProfile::collect(const NoiseProfile& o) {
	for (int i = 0; i < NCODES; i++)
		p[i] += o.p[i];
}

void NoiseProfile::read(FILE *fi) {
	int tmp_ncodes;

	memset(c, 0, sizeof(c));
	assert(fscanf(fi, "%d", &tmp_ncodes) == 1);
	assert(tmp_ncodes == NCODES);
	for (int i = 0; i < NCODES; i++)
	  assert(fscanf(fi, "%lf", &p[i]) == 1);
}

void NoiseProfile::write(FILE *fo) {
	fprintf(fo, "%d\n", NCODES);
	for (int i = 0; i < NCODES - 1; i++) {
		fprintf(fo, "%.10g ", p[i]);
	}
	fprintf(fo, "%.10g\n", p[NCODES - 1]);
}

void NoiseProfile::startSimulation() {
	pc = new double[NCODES];

	for (int i = 0; i < NCODES; i++) {
		pc[i] = p[i];
		if (i > 0) pc[i] += pc[i - 1];
	}
}

std::string NoiseProfile::simulate(simul* sampler, int len) {
	std::string readseq = "";

	for (int i = 0; i < len; i++) {
		readseq.push_back(getCharacter(sampler->sample(pc, NCODES)));
	}
	return readseq;
}

void NoiseProfile::finishSimulation() {
	delete[] pc;
}

#endif /* NOISEPROFILE_H_ */
