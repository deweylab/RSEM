/* Copyright (c) 2017
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#include <new>
#include <cmath>
#include <limits>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "QualDist.hpp"

QualDist::QualDist(int mode) : mode(mode), p_init(NULL), ss_init(NULL), p_tran(NULL), ss_tran(NULL) {
	p_init = new double[QSIZE];
	p_tran = new double[QSIZE][QSIZE];

	if (mode == 0) {
		ss_init = new double[QSIZE];
		ss_tran = new double[QSIZE][QSIZE];		
	}
}

QualDist::~QualDist() {
	if (p_init != NULL) delete[] p_init;
	if (p_tran != NULL) delete[] p_tran;
	if (ss_init != NULL) delete[] ss_init;
	if (ss_tran != NULL) delete[] ss_tran;
}

void QualDist::clear() {
	memset(p_init, 0, sizeof(double) * QSIZE);
	memset(p_tran, 0, sizeof(double) * QSIZE * QSIZE);
}

void QualDist::collect(const QualDist* o) {
	for (int i = 0; i < QSIZE; ++i) p_init[i] += o->p_init[i];
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < QSIZE; ++j) p_tran[i][j] += o->p_tran[i][j];
}

void QualDist::finish() {
	memcpy(ss_init, p_init, sizeof(double) * QSIZE);
	memcpy(ss_tran, p_tran, sizeof(double) * QSIZE * QSIZE);
	ss2p();
	p2logp();
}

double QualDist::calcLogP() const {
	double logp = 0.0;

	for (int i = 0; i < QSIZE; ++i)
		if (ss_init[i] > 0.0) logp += ss_init[i] * p_init[i];

	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < QSIZE; ++j)
			if (ss_tran[i][j] > 0.0) logp += ss_tran[i][j] * p_tran[i][j];

	return logp;
}

void QualDist::read(std::ifstream& fin, int choice) {
	std::string line;
	int tmp_qsize;
	double *in_init = NULL, (*in_tran)[QSIZE] = NULL;

	switch(choice) {
		case 0: in_init = p_init; in_tran = p_tran; break;
		case 1: in_init = ss_init; in_tran = ss_tran;
	}

	assert((fin>> tmp_qsize) && (tmp_qsize == QSIZE));
	for (int i = 0; i < QSIZE; ++i) assert(fin>> in_init[i]);
	for (int i = 0; i < QSIZE; ++i) 
		for (int j = 0; j < QSIZE; ++j) assert(fin>> in_tran[i][j]);
	getline(fin, line);

	if (mode == 0 && choice == 0) p2logp();
	if (mode == 2) prepare_for_simulation();
}

void QualDist::write(std::ofstream& fout, int choice) {
	double *out_init = NULL, (*out_tran)[QSIZE] = NULL;

	switch(choice) {
		case 0: ss2p(); out_init = p_init; out_tran = p_tran; break;
		case 1: out_init = ss_init; out_tran = ss_tran;
	}

	fout<< "#qd\tQualDist\t"<< QSIZE + 5<< "\tformat: SIZE; P_init, one line with SIZE values; P_tran, a SIZE x SIZE matrix"<< std::endl;
	fout<< QSIZE<< std::endl;  
	for (int i = 0; i < QSIZE - 1; ++i) fout<< out_init[i]<< '\t';
	fout<< out_init[QSIZE - 1]<< std::endl<< std::endl;
	for (int i = 0; i < QSIZE; ++i) {
		for (int j = 0; j < QSIZE - 1 ; ++j) fout<< out_tran[i][j]<< '\t';
		fout<< out_tran[i][QSIZE - 1]<< std::endl;
	}
	fout<< std::endl<< std::endl;
}

void QualDist::ss2p() {
	double sum = 0.0;
	for (int i = 0; i < QSIZE; ++i) sum += ss_init[i];
	assert(sum > 0.0);
	for (int i = 0; i < QSIZE; ++i) p_init[i] = ss_init[i] / sum;

	for (int i = 0; i < QSIZE; ++i) {
		sum = 0.0;
		for (int j = 0; j < QSIZE; ++j) sum += ss_tran[i][j];
		if (sum > 0.0)
			for (int j = 0; j < QSIZE; ++j) p_tran[i][j] = ss_tran[i][j] / sum;
		else memset(p_tran[i], 0, sizeof(double) * QSIZE);
	}
}

void QualDist::p2logp() {
	for (int i = 0; i < QSIZE; ++i) p_init[i] = (p_init[i] > 0.0 ? log(p_init[i]) : -std::numeric_limits<double>::infinity());
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < QSIZE; ++j)
			p_tran[i][j] = (p_tran[i][j] > 0.0 ? log(p_tran[i][j]) : -std::numeric_limits<double>::infinity());	
}

void QualDist::prepare_for_simulation() {
	for (int i = 1; i < QSIZE; ++i) p_init[i] += p_init[i - 1];
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 1; j < QSIZE; ++j) 
			p_tran[i][j] += p_tran[i][j - 1];
}
