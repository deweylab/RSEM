#ifndef QUALDIST_H_
#define QUALDIST_H_

#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>

#include "simul.h"

//from 33 to 126 to encode 0 to 93
class QualDist {
public:
	QualDist() {
		memset(p_init, 0, sizeof(p_init));
		memset(p_tran, 0, sizeof(p_tran));
	}

	QualDist& operator=(const QualDist&);

	void update(const std::string&);
	void finish();

	double getProb(const std::string&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int);
	void finishSimulation();

private:
	static const int SIZE = 100;

	double p_init[SIZE];
	double p_tran[SIZE][SIZE]; //p_tran[a][b] = p(b|a)

	int c2q(char c) { assert(c >= 33 && c <= 126); return c - 33; }

	double *qc_init, (*qc_trans)[SIZE];
	char q2c(int qval) { return (char)(qval + 33); }
};

QualDist& QualDist::operator=(const QualDist& rv) {
	if (this == &rv) return *this;

	memcpy(p_init, rv.p_init, sizeof(rv.p_init));
	memcpy(p_tran, rv.p_tran, sizeof(rv.p_tran));

	return *this;
}

void QualDist::update(const std::string& qual) {
	int len = qual.size();

	assert(len > 0);
	++p_init[c2q(qual[0])];

	for (int i = 1; i < len; i++) {
		++p_tran[c2q(qual[i - 1])][c2q(qual[i])];
	}
}

void QualDist::finish() {
	double sum;

	sum = 0.0;
	for (int i = 0; i < SIZE; i++) sum += p_init[i];
	for (int i = 0; i < SIZE; i++) p_init[i] /= sum;

	for (int i = 0; i < SIZE; i++) {
		sum = 0.0;
		for (int j = 0; j < SIZE; j++) sum += p_tran[i][j];
		if (sum <= 0.0) continue;
		//if (isZero(sum)) continue;
		for (int j = 0; j < SIZE; j++) p_tran[i][j] /= sum;
	}
}

double QualDist::getProb(const std::string& qual) {
	int len = qual.size();
	double prob = 1.0;

	assert(len > 0);
	prob *= p_init[c2q(qual[0])];
	for (int i = 1; i < len; i++) {
		prob *= p_tran[c2q(qual[i - 1])][c2q(qual[i])];
	}

	return prob;
}

void QualDist::read(FILE *fi) {
	int tmp_size;

	assert(fscanf(fi, "%d", &tmp_size) == 1);
	assert(tmp_size == SIZE);

	for (int i = 0; i < SIZE; i++) { assert(fscanf(fi, "%lf", &p_init[i]) == 1); }
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) { assert(fscanf(fi, "%lf", &p_tran[i][j]) == 1); }
	}
}

void QualDist::write(FILE *fo) {
	fprintf(fo, "%d\n", SIZE);
	for (int i = 0; i < SIZE - 1; i++) { fprintf(fo, "%.10g ", p_init[i]); }
	fprintf(fo, "%.10g\n", p_init[SIZE - 1]);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE -1 ; j++) fprintf(fo, "%.10g ", p_tran[i][j]);
		fprintf(fo, "%.10g\n", p_tran[i][SIZE - 1]);
	}
}

void QualDist::startSimulation() {
	qc_init = new double[SIZE];
	qc_trans = new double[SIZE][SIZE];

	for (int i = 0; i < SIZE; i++) {
		qc_init[i] = p_init[i];
		if (i > 0) qc_init[i] += qc_init[i - 1];
	}

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++) {
			qc_trans[i][j] = p_tran[i][j];
			if (j > 0) qc_trans[i][j] += qc_trans[i][j - 1];
		}
}

std::string QualDist::simulate(simul* sampler, int len) {
	int qval, old_qval;
	std::string qual = "";

	qval = sampler->sample(qc_init, SIZE);
	qual.push_back(q2c(qval));
	for (int i = 1; i < len; i++) {
		old_qval = qval;
		qval = sampler->sample(qc_trans[old_qval], SIZE);
		qual.push_back(q2c(qval));
	}

	return qual;
}

void QualDist::finishSimulation() {
	delete[] qc_init;
	delete[] qc_trans;
}
#endif /* QUALDIST_H_ */
