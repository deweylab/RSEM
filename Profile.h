#ifndef PROFILE_H_
#define PROFILE_H_

#include<cstdio>
#include<cstring>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "simul.h"


class Profile {
public:
	Profile(int = 1000);
	~Profile() {
		delete[] p;
	}

	Profile& operator=(const Profile&);

	void init();
	void update(const std::string&, const RefSeq&, int, int, double);
	void finish();

	double getProb(const std::string&, const RefSeq&, int, int);

	void collect(const Profile&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int, int, int, const RefSeq&);
	void finishSimulation();

private:
	static const int NCODES = 5;

	int proLen; // profile length
	int size; // # of items in p;
	double (*p)[NCODES][NCODES]; //profile matrices

	double (*pc)[NCODES][NCODES]; // for simulation
};

Profile::Profile(int maxL) {
	proLen = maxL;
	size = proLen * NCODES * NCODES;
	p = new double[proLen][NCODES][NCODES];
	memset(p, 0, sizeof(double) * size);

	//set initial parameters
	int N = NCODES - 1;
	double probN = 1e-5, portionC = 0.99; //portionC, among ACGT, the portion of probability mass the correct base takes
	double probC, probO;

	for (int i = 0; i < proLen; i++) {
		for (int j = 0; j < NCODES - 1; j++) {
			p[i][j][N] = probN;
			probC = portionC * (1.0 - probN);
			probO = (1.0 - portionC) / (NCODES - 2) * (1.0 - probN);

			for (int k = 0; k < NCODES - 1; k++) {
				p[i][j][k] = (j == k ? probC : probO);
			}
		}
		p[i][N][N] = probN;
		for (int k = 0; k < NCODES - 1; k++)
				p[i][N][k] = (1.0 - probN) / (NCODES - 1);
	}
}

Profile& Profile::operator=(const Profile& rv) {
	if (this == &rv) return *this;
	if (proLen != rv.proLen) {
		delete[] p;
		proLen = rv.proLen;
		size = rv.size;
		p = new double[rv.proLen][NCODES][NCODES];
	}
	memcpy(p, rv.p, sizeof(double) * rv.size);

	return *this;
}

void Profile::init() {
	memset(p, 0, sizeof(double) * size);
}

void Profile::update(const std::string& readseq, const RefSeq& refseq, int pos, int dir, double frac) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		p[i][refseq.get_id(i + pos, dir)][get_base_id(readseq[i])] += frac;
	}
}

void Profile::finish() {
	double sum;

	for (int i = 0; i < proLen; i++) {
		for (int j = 0; j < NCODES; j++) {
			sum = 0.0;
			for (int k = 0; k < NCODES; k++) sum += p[i][j][k];
			if (sum < EPSILON) {
				for (int k = 0; k < NCODES; k++) p[i][j][k] = 0.0;
				continue;
			}
			for (int k = 0; k < NCODES; k++) p[i][j][k] /= sum;
		}
	}
}

double Profile::getProb(const std::string& readseq, const RefSeq& refseq, int pos, int dir) {
	double prob = 1.0;
	int len = readseq.size();

	for (int i = 0; i < len; i++) {
		prob *= p[i][refseq.get_id(i + pos, dir)][get_base_id(readseq[i])];

	}

	return prob;
}

void Profile::collect(const Profile& o) {
	for (int i = 0; i < proLen; i++)
		for (int j = 0; j < NCODES; j++)
			for (int k = 0; k < NCODES; k++)
				p[i][j][k] += o.p[i][j][k];
}

void Profile::read(FILE *fi) {
	int tmp_prolen, tmp_ncodes;
	assert(fscanf(fi, "%d %d", &tmp_prolen, &tmp_ncodes) == 2);
	assert(tmp_ncodes == NCODES);
	if (tmp_prolen != proLen) {
		delete[] p;
		proLen = tmp_prolen;
		size = proLen * NCODES * NCODES;
		p = new double[proLen][NCODES][NCODES];
		memset(p, 0, sizeof(double) * size);
	}

	for (int i = 0; i < proLen; i++)
		for (int j = 0; j < NCODES; j++)
			for (int k = 0; k < NCODES; k++)
			  assert(fscanf(fi, "%lf", &p[i][j][k]) == 1);
}

void Profile::write(FILE* fo) {
	fprintf(fo, "%d %d\n", proLen, NCODES);
	for (int i = 0; i < proLen; i++) {
		for (int j = 0; j < NCODES; j++) {
			for (int k = 0; k < NCODES - 1; k++)
				fprintf(fo, "%.10g ", p[i][j][k]);
			fprintf(fo, "%.10g\n", p[i][j][NCODES - 1]);
		}
		if (i < proLen - 1) { fprintf(fo, "\n"); }
	}
}

void Profile::startSimulation() {
	pc = new double[proLen][NCODES][NCODES];
	for (int i = 0; i < proLen; i++) {
		for (int j = 0; j < NCODES; j++)
			for (int k = 0; k < NCODES; k++) {
				pc[i][j][k] = p[i][j][k];
				if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
			}
		//avoid sampling from 0.0!!!
		double cp_sum, cp_d, cp_n;
		cp_sum = cp_d = cp_n = 0.0;

		for (int j = 0; j < NCODES - 1; j++) {
		  cp_sum += pc[i][j][NCODES - 1];
		  cp_d += p[i][j][j];
		  cp_n += p[i][j][NCODES - 1];
		}

		if (cp_sum == 0.0) continue;

		double p_d, p_o, p_n;
		p_d = cp_d / cp_sum;
		p_n = cp_n / cp_sum;
		p_o = (1.0 - p_d - p_n) / (NCODES - 2);
		for (int j = 0; j < NCODES - 1; j++) {
		  if (pc[i][j][NCODES - 1] > 0.0) continue;
		  for (int k = 0; k < NCODES; k++) {
		    if (k == j) pc[i][j][k] = p_d;
		    else if (k == NCODES - 1) pc[i][j][k] = p_n;
		    else pc[i][j][k] = p_o;
		    if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
		  }
		}
		if (pc[i][NCODES - 1][NCODES - 1] == 0.0) {
		  p_o = (1.0 - p_n) / (NCODES - 1);
		  for (int k = 0; k < NCODES; k++) {
		    pc[i][NCODES - 1][k] = (k < NCODES - 1 ? p_o : p_n);
		    if (k > 0) pc[i][NCODES - 1][k] += pc[i][NCODES - 1][k - 1];
		  }
		}
	}

}

std::string Profile::simulate(simul* sampler, int len, int pos, int dir, const RefSeq& refseq) {
	std::string readseq = "";

	for (int i = 0; i < len; i++) {
		readseq.push_back(getCharacter(sampler->sample(pc[i][refseq.get_id(i + pos, dir)], NCODES)));
	}
	return readseq;
}

void Profile::finishSimulation() {
	delete[] pc;
}

#endif /* PROFILE_H_ */
