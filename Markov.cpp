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
#include <vector>
#include <string>
#include <fstream>

#include "utils.h"
#include "Markov.hpp"

const char Markov::state2char[Markov::NSTATES] = {0, 'S', 'M', 'I', 'D', 'S'};

const int Markov::lens[Markov::NSTATES] = {2, 2, 4, 3, 3, 1};
const int Markov::offset[Markov::NSTATES] = {1, 1, 2, 2, 2, 5};

const int Markov::s2pos[Markov::NSTATES] = {-1, 0, -1, 1, -1, 2};

Markov::Markov(bool to_log_space) : to_log_space(to_log_space) {
	memset(P, 0, sizeof(P));
	memset(probB, 0, sizeof(probB));
	ss_P = NULL; ss_probB = NULL;
}

Markov::~Markov() {
	if (ss_P != NULL) delete[] ss_P;
	if (ss_probB != NULL) delete[] ss_probB;
}

void Markov::collect(const Markov* o) {
	for (int i = 0; i < NSTATES; ++i)
		for (int j = 0; j < NSTATES; ++j)
			P[i][j] += o->P[i][j];

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < NCODES; ++j)
			probB[i][j] += o->probB[i][j];
}

void Markov::finish() {
	// collect sufficient statistics
	if (ss_P == NULL) ss_P = new double[NSTATES][NSTATES];
	memcpy(ss_P, P, sizeof(P));
	if (ss_probB == NULL) ss_probB = new double[3][NCODES];
	memcpy(ss_probB, probB, sizeof(probB));
	// calculate probabilities
  	ss2p();
	// convert to log space
  	if (to_log_space) {
		for (int i = 0; i < NSTATES; ++i)
			for (int j = 0; j < NSTATES; ++j)
				P[i][j] = (P[i][j] == 0.0 ? -1000.0 : log(P[i][j]));
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < NCODES; ++j)
				probB[i][j] = (probB[i][j] == 0.0 ? -1000.0 : log(probB[i][j]));
	}
}

void Markov::clear() {
	memset(P, 0, sizeof(P));
	memset(probB, 0, sizeof(probB));
}

void Markov::read(std::ifstream& fin) {
	std::string line;
	int tmp_nstates, tmp_ncodes;

	assert((fin>> tmp_nstates>> tmp_ncodes) && (tmp_nstates == NSTATES) && (tmp_ncodes == NCODES));

	for (int i = 0; i < NSTATES; ++i)
		for (int j = 0; j < NSTATES; ++j)
			assert(fin>> P[i][j]);

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < NCODES; ++j)
			assert(fin>> probB[i][j]);

	getline(fin, line);
}

void Markov::write(std::ofstream& fout, bool isProb) {
	if (to_log_space) ss2p();

	fout<< "#markov "<< 1 + (NSTATES + 1) + 5<< " Markov model for CIGARs, format: NSTATES NCODES; P; probB"<< std::endl;

	fout<< NSTATES<< '\t'<< NCODES<< std::endl;
  
	for (int i = 0; i < NSTATES; ++i) {
		for (int j = 0; j < NSTATES - 1; ++j) fout<< (isProb ? P[i][j] : ss_P[i][j])<< '\t';
			fout<< (isProb ? P[i][NSTATES - 1] : ss_P[i][NSTATES － 1])<< std::endl;
	}
	fout<< std::endl;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < NCODES - 1; ++j) fout<< (isProb ? probB[i][j] : ss_probB[i][j])<< '\t';
		fout<< (isProb ? probB[i][NCODES - 1] : ss_probB[i][NCODES - 1])<< std::endl;
	}
	fout<< std::endl<< std::endl;
}

void Markov::prepare_for_simulation() {
	if (P[M][S2] > 0.0 && P[S2][S2] == 0.0) P[S2][S2] = 1.0; // in case S2 only appear as the last base of read/mate

	for (int i = 0; i < NSTATES; ++i) {
		P[i][0] = P[i][offset[i]];
		for (int j = 1; j < lens[j]; ++j)
			P[i][j] = P[i][j - 1] + P[i][offset[i] + j];
	}

	for (int i = 0; i < 3; ++i)
		for (int j = 1; j < NCODES; ++j)
			probB[i][j] += probB[i][j - 1];
}

void Markov::ss2p() {
	double sum;

	if (!to_log_space) {
		for (int i = 0; i < NSTATES; ++i)
			for (int j = 0; j < NSTATES; ++j)
				ss_P[i][j] += 0.1; // add a small pseudo count for pre-estimation with unique alignments only
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < NCODES; ++j)
				ss_probB[i][j] += pseudoC[j];
	}

	for (int i = 0; i < NSTATES; ++i) {
		sum = 0.0;
		for (int j = 0; j < NSTATES; ++j) sum += ss_P[i][j];
		if (sum > 0.0)
			for (int j = 0; j < NSTATES; ++j) P[i][j] = ss_P[i][j] / sum;
		else
			memset(P[i], 0, sizeof(double) * NSTATES);
	}

	for (int i = 0; i < 3; ++i) {
		sum = 0.0;
		for (int j = 0; j < NCODES; ++j) sum += ss_probB[i][j];
		if (sum > 0.0)
			for (int j = 0; j < NCODES; ++j) probB[i][j] = ss_probB[i][j] / sum;
		else
			memset(probB[i], 0, sizeof(double) * NCODES);
	}
}
