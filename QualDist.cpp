/* Copyright (c) 2016
   Bo Li (University of California, Berkeley)
   bli25@berkeley.edu

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

#include <cmath>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "QualDist.hpp"

QualDist::QualDist() {
	memset(p_init, 0, sizeof(p_init));
	memset(p_tran, 0, sizeof(p_tran));

	ss_init = NULL; ss_tran = NULL;
}

QualDist::~QualDist() {
	if (ss_init != NULL) delete[] ss_init;
	if (ss_tran != NULL) delete[] ss_tran;
}

void QualDist::collect(const QualDist* o) {
	for (int i = 0; i < QSIZE; ++i) p_init[i] += o->p_init[i];
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < QSIZE; ++j) p_tran[i][j] += o->p_tran[i][j];
}

void QualDist::finish() {
	double sum, value;

	// collect sufficient statistics
	if (ss_init == NULL) ss_init = new double[QSIZE];
	memcpy(ss_init, p_init, sizeof(double) * QSIZE);
	if (ss_trans == NULL) ss_tran = new double[QSIZE][QSIZE];
	memcpy(ss_tran, p_tran, sizeof(double) * QSIZE * QSIZE);

	// calculate QualDist
	sum = 0.0;
	for (int i = 0; i < QSIZE; ++i) sum += p_init[i];
	assert(sum > 0.0);
	for (int i = 0; i < QSIZE; ++i) p_init[i] /= sum;

	for (int i = 0; i < QSIZE; ++i) {
		sum = 0.0;
		for (int j = 0; j < QSIZE; ++j) sum += p_tran[i][j];
		if (sum > 0.0)
			for (int j = 0; j < QSIZE; ++j) p_tran[i][j] /= sum;
	}
}

double QualDist::calcLogP() const {
	double logp = 0.0;

	for (int i = 0; i < QSIZE; ++i)
		if (ss_init[i] > 0.0) logp += ss_init[i] * log(p_init[i]);

	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < QSIZE; ++j)
			if (ss_tran[i][j] > 0.0) logp += ss_tran[i][j] * log(p_tran[i][j]);

	return logp;
}

void QualDist::read(std::ifstream& fin) {
	std::string line;
	int tmp_qsize;

	assert((fin>> tmp_qsize) && (tmp_qsize == QSIZE));
	for (int i = 0; i < QSIZE; ++i) assert(fin>> p_init[i]);
	for (int i = 0; i < QSIZE; ++i) 
		for (int j = 0; j < QSIZE; ++j) assert(fin>> p_tran[i][j]);

	getline(fin, line);
}

void QualDist::write(std::ofstream& fout, bool isProb) {
	// comment sign, model key word, number of lines, description and format
	fout<< "#qd "<< QSIZE + 5<< " QualDist, format: SIZE; P_init; P_tran, SIZE lines"<< std::endl;

	fout<< QSIZE<< std::endl;  
	for (int i = 0; i < QSIZE - 1; ++i) fout<< (isProb ? p_init[i] : ss_init[i])<< '\t';
	fout<< (isProb ? p_init[QSIZE - 1] : ss_init[QSIZE - 1])<< std::endl<< std::endl;
	for (int i = 0; i < QSIZE; ++i) {
		for (int j = 0; j < QSIZE - 1 ; ++j) fout<< (isProb ? p_tran[i][j] : ss_tran[i][j])<< '\t';
		fout<< (isProb ? p_tran[i][QSIZE - 1] : ss_tran[i][QSIZE - 1])<< std::endl;
	}
	fout<< std::endl<< std::endl;
}

void QualDist::prepare_for_simulation() {
	for (int i = 1; i < QSIZE; ++i) p_init[i] += p_init[i - 1];

	for (int i = 0; i < QSIZE; ++i)
		for (int j = 1; j < QSIZE; ++j) 
			p_tran[i][j] += p_tran[i][j - 1];
}
