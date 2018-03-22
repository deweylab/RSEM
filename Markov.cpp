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
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "Markov.hpp"

const char Markov::state2char[Markov::NSTATES] = {0, 'S', 'M', 'I', 'D', 'S'};

const int Markov::lens[Markov::NSTATES] = {2, 2, 4, 3, 3, 1};
const int Markov::offset[Markov::NSTATES] = {1, 1, 2, 2, 2, 5};

const int Markov::s2pos[Markov::NSTATES] = {-1, 0, -1, 1, -1, 2};

Markov::Markov(model_mode_type mode) : mode(mode), P(NULL), ss_P(NULL), probB(NULL), ss_probB(NULL) {
	P = new double[NSTATES][NSTATES];
	probB = new double[3][NCODES];
	if (mode == FIRST_PASS || mode == MASTER) {
		ss_P = new double[NSTATES][NSTATES];
		ss_probB = new double[3][NCODES];
	}
	if (mode == INIT) init();
}

Markov::~Markov() {
	if (P != NULL) delete[] P;
	if (probB != NULL) delete[] probB;
	if (ss_P != NULL) delete[] ss_P;
	if (ss_probB != NULL) delete[] ss_probB;
}

void Markov::clear() {
	memset(P, 0, sizeof(double) * NSTATES * NSTATES);
	memset(probB, 0, sizeof(double) * 3 * NCODES);
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
	memcpy(ss_P, P, sizeof(double) * NSTATES * NSTATES);
	memcpy(ss_probB, probB, sizeof(double) * 3 * NCODES);
  	ss2p();
  	if (mode == MASTER) p2logp();
}

void Markov::read(std::ifstream& fin, int choice) {
	std::string line;
	int tmp_nstates, tmp_ncodes;
	double (*in_P)[NSTATES] = NULL, (*in_probB)[NCODES] = NULL;

	switch(choice) {
		case 0: in_P = P; in_probB = probB; break;
		case 1: in_P = ss_P; in_probB = ss_probB;
	}

	assert((fin>> line) && (line == "#markov"));
	assert(getline(fin, line));
	assert((fin>> tmp_nstates>> tmp_ncodes) && (tmp_nstates == NSTATES) && (tmp_ncodes == NCODES));
	for (int i = 0; i < NSTATES; ++i)
		for (int j = 0; j < NSTATES; ++j)
			assert(fin>> in_P[i][j]);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < NCODES; ++j)
			assert(fin>> in_probB[i][j]);
	assert(getline(fin, line));

	if (mode == MASTER && choice == 0) p2logp();
	if (mode == SIMULATION) prepare_for_simulation();
}

void Markov::write(std::ofstream& fout, int choice) {
	double (*out_P)[NSTATES] = NULL, (*out_probB)[NCODES] = NULL;

	switch(choice) {
		case 0: if (mode == MASTER) ss2p(); out_P = P; out_probB = probB; break;
		case 1: out_P = ss_P; out_probB = ss_probB;
	}

	fout<< "#markov\tMarkov model for CIGARs\t"<< 1 + (NSTATES + 1) + 5<< "\tformat: NSTATES NCODES; P, NSTATES x NSTATES; probB, 3 x NCODES"<< std::endl;
	fout<< NSTATES<< '\t'<< NCODES<< std::endl;
	for (int i = 0; i < NSTATES; ++i) {
		for (int j = 0; j < NSTATES - 1; ++j) fout<< out_P[i][j]<< '\t';
			fout<< out_P[i][NSTATES - 1]<< std::endl;
	}
	fout<< std::endl;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < NCODES - 1; ++j) fout<< out_probB[i][j]<< '\t';
		fout<< out_probB[i][NCODES - 1]<< std::endl;
	}
	fout<< std::endl<< std::endl;
}

void Markov::init() {
	memset(P, 0, sizeof(double) * NSTATES * NSTATES);
	P[Begin][0] = 0.9; P[Begin][1] = 0.1;
	P[S1][0] = 0.1; P[S1][1] = 0.9;
	P[M][0] = 0.7; P[M][1] = P[M][2] = P[M][3] = 0.1;
	P[I][0] = 0.8; P[I][1] = 0.15; P[I][2] = 0.05;
	P[D][0] = 0.8; P[D][1] = 0.05; P[D][2] = 0.15;
	P[S2][0] = 1.0;
	for (int i = 0; i < 3; ++i) 
		for (int j = 0; j < NCODES; ++j)
			probB[i][j] = prior_noise[j];
	p2logp();
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

void Markov::p2logp() {
	for (int i = 0; i < NSTATES; ++i)
		for (int j = 0; j < NSTATES; ++j)
				P[i][j] = (P[i][j] > 0.0 ? log(P[i][j]) : NEGINF);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < NCODES; ++j)
			probB[i][j] = (probB[i][j] > 0.0 ? log(probB[i][j]) : NEGINF);
}
