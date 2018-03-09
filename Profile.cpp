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
#include "Profile.hpp"

Profile::Profile(model_mode_type mode, int maxL) : mode(mode), maxL(maxL), p(NULL), ss(NULL) {
	assert(maxL > 0);
	p = new double[maxL][NCODES][NCODES];
	if (mode == FIRST_PASS || mode == MASTER) ss = new double[maxL][NCODES][NCODES];
	if (mode == INIT) init();
}

Profile::~Profile() {
	if (p != NULL) delete[] p;
	if (ss != NULL) delete[] ss;
}

void Profile::clear() {
	memset(p, 0, sizeof(double) * maxL * NCODES * NCODES);
}

void Profile::collect(const Profile* o) {
	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] += o->p[i][j][k];
}

void Profile::finish() {
	memcpy(ss, p, sizeof(double) * maxL * NCODES * NCODES);
	ss2p(); // calculate p
	if (mode == MASTER) p2logp(); // convert to log space
}

void Profile::read(std::ifstream& fin, int choice) {
	std::string line;
	int tmp_maxl, tmp_ncodes;
	double (*in)[NCODES][NCODES] = NULL;

	switch(choice) {
		case 0: in = p; break;
		case 1: in = ss;
	}

	assert((fin>> line) && (line == "#pro"));
	assert(getline(fin, line));
	assert((fin>> tmp_maxl>> tmp_ncodes) && (tmp_maxl == maxL) && (tmp_ncodes == NCODES));
	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				assert(fin>> in[i][j][k]);
	assert(getline(fin, line));

	if (mode == MASTER && choice == 0) p2logp();
	if (mode == SIMULATION) prepare_for_simulation();
}

void Profile::write(std::ofstream& fout, int choice) {
	double (*out)[NCODES][NCODES] = NULL;

	switch(choice) {
		case 0: if (mode == MASTER) ss2p(); out = p; break;
		case 1: out = ss;
	}

	fout<< "#pro\tProfile\t"<< 1 + (NCODES + 1) * maxL + 1<< "\tformat: maxL NCODES; P[POS][REF_BASE][OBSERVED_BASE], maxL blocks separated by a blank line, each block is a NCODES x NCODES matrix"<< std::endl;
	fout<< maxL<< '\t'<< NCODES<< std::endl;
	for (int i = 0; i < maxL; ++i) {
		for (int j = 0; j < NCODES; ++j) {
			for (int k = 0; k < NCODES - 1; ++k)
				fout<< out[i][j][k]<< '\t';
			fout<< out[i][j][NCODES - 1]<< std::endl;
		}
		fout<< std::endl;
	}
	fout<< std::endl;
}

void Profile::init() {
	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] = (j < NCODES - 1 ? (j == k ? prior_aligned[0] : (k < NCODES - 1 ? prior_aligned[1] : prior_aligned[2])): prior_aligned[3]);
	p2logp();
}

void Profile::prepare_for_simulation() {
	for (int i = 0; i < maxL; ++i) 
		for (int j = 0; j < NCODES; ++j)
			for (int k = 1; k < NCODES; ++k)
				p[i][j][k] += p[i][j][k - 1];
}

void Profile::ss2p() {
	double sum;	
	for (int i = 0; i < maxL; ++i) 
		for (int j = 0; j < NCODES; ++j) {
			sum = 0.0;
			for (int k = 0; k < NCODES; ++k) {
				p[i][j][k] = ss[i][j][k] + pseudo_count * (j < NCODES - 1 ? (j == k ? prior_aligned[0] : (k < NCODES - 1 ? prior_aligned[1] : prior_aligned[2])): prior_aligned[3]);
				sum += p[i][j][k];
			}
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] /= sum;
		}
}

void Profile::p2logp() {
	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] = log(p[i][j][k]);
}
