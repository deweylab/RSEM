#ifndef QPROFILE_H_
#define QPROFILE_H_

#include<cstdio>
#include<cstring>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "simul.h"


class QProfile {
public:
	QProfile();
	QProfile& operator=(const QProfile&);

	void init();
	void update(const std::string&, const std::string&, const RefSeq&, int, int, double);
	void finish();

	double getProb(const std::string&, const std::string&, const RefSeq&, int, int);

	void collect(const QProfile&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int, int, int, const std::string&, const RefSeq&);
	void finishSimulation();

private:
	static const int NCODES = 5; // number of possible codes
	static const int SIZE = 100;

	double p[SIZE][NCODES][NCODES]; // p[q][r][c] = p(c|r,q)

	//make sure that quality score in [0, 93]
	int c2q(char c) { assert(c >= 33 && c <= 126); return c - 33; }

	double (*pc)[NCODES][NCODES]; // for simulation
};

QProfile::QProfile() {
	memset(p, 0, sizeof(p));

	//make initialized parameters
	//ASSUME order of A, C, G, T, N
	int N = NCODES - 1;
	double probN = 1e-5;
	double probC, probO; // current, other

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < NCODES - 1; j++) {
			p[i][j][N] = probN;

			probO = exp(-i / 10.0 * log(10.0));
			probC = 1.0 - probO;
			probO /= (NCODES - 2);

			probC *= (1.0 - probN);
			probO *= (1.0 - probN);

			assert(probC >= 0.0 && probO >= 0.0);

			for (int k = 0; k < NCODES - 1; k++) {
				if (j == k) p[i][j][k] = probC;
				else p[i][j][k] = probO;
			}
		}
		p[i][N][N] = probN;
		for (int k = 0; k < NCODES - 1; k++)
			p[i][N][k] = (1.0 - probN) / (NCODES - 1);
	}
}

QProfile& QProfile::operator=(const QProfile& rv) {
	if (this == &rv) return *this;
	memcpy(p, rv.p, sizeof(rv.p));
	return *this;
}

void QProfile::init() {
	memset(p, 0, sizeof(p));
}

void QProfile::update(const std::string& readseq, const std::string& qual, const RefSeq& refseq, int pos, int dir, double frac) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
	  p[c2q(qual[i])][refseq.get_id(i + pos, dir)][get_base_id(readseq[i])] += frac;
	}
}

void QProfile::finish() {
	double sum;

	for (int i = 0; i < SIZE; i++) {
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

double QProfile::getProb(const std::string& readseq, const std::string& qual, const RefSeq& refseq, int pos, int dir) {
	double prob = 1.0;
	int len = readseq.size();

	for (int i = 0; i < len; i++) {
		prob *= p[c2q(qual[i])][refseq.get_id(i + pos, dir)][get_base_id(readseq[i])];
	}

	return prob;
}

void QProfile::collect(const QProfile& o) {
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < NCODES; j++)
			for (int k = 0; k < NCODES; k++)
				p[i][j][k] += o.p[i][j][k];
}

void QProfile::read(FILE *fi) {
	int tmp_size, tmp_ncodes;
	assert(fscanf(fi, "%d %d", &tmp_size, &tmp_ncodes) == 2);
	assert(tmp_size == SIZE && tmp_ncodes == NCODES);
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < NCODES; j++)
			for (int k = 0; k < NCODES; k++)
			  assert(fscanf(fi, "%lf", &p[i][j][k]) == 1);
}

void QProfile::write(FILE *fo) {
	fprintf(fo, "%d %d\n", SIZE, NCODES);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < NCODES; j++) {
			for (int k = 0; k < NCODES - 1; k++)
				fprintf(fo, "%.10g ", p[i][j][k]);
			fprintf(fo, "%.10g\n", p[i][j][NCODES - 1]);
		}
		if (i < SIZE - 1) { fprintf(fo, "\n"); }
	}
}

void QProfile::startSimulation() {
	pc = new double[SIZE][NCODES][NCODES];
	for (int i = 0; i < SIZE; i++) {
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

std::string QProfile::simulate(simul* sampler, int len, int pos, int dir, const std::string& qual, const RefSeq& refseq) {
	std::string readseq = "";

	for (int i = 0; i < len; i++) {
		readseq.push_back(getCharacter(sampler->sample(pc[c2q(qual[i])][refseq.get_id(i + pos, dir)], NCODES)));
	}
	return readseq;
}

void QProfile::finishSimulation() {
	delete[] pc;
}

#endif /* QPROFILE_H_ */
