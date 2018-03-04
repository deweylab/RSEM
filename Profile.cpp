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
#include <cstdlib>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "Profile.hpp"

Profile::Profile(bool to_log_space, int maxL) : to_log_space(to_log_space), maxL(maxL) {
	p = NULL; ss = NULL;
	if (maxL > 0) allocate_space(maxL);
}

Profile::~Profile() {
	if (p != NULL) { 
		for (int i = 0; i < maxL; ++i) free(p[i]);
		free(p); 
	}
	if (ss != NULL) delete[] ss;
}

void Profile::collect(const Profile* o) {
	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] += o->p[i][j][k];
}

void Profile::finish(int length) {
	// expand p if necessary, length is the maximum mate length observed
	expand(length);
	// collect sufficient statistics
	if (ss == NULL) ss = new double[maxL][NCODES][NCODES];
	memcpy(ss, p, sizeof(double) * maxL * NCODES * NCODES);
	// calculate p
	ss2p();
	// convert p into log space
	if (to_log_space) 
		for (int i = 0; i < maxL; ++i)
			for (int j = 0; j < NCODES; ++j)
				for (int k = 0; k < NCODES; ++k)
					p[i][j][k] = log(p[i][j][k]);
}

void Profile::clear() {
	size_t num = sizeof(double) * NCODES * NCODES;
	for (int i = 0; i < maxL; ++i) memset(p[i], 0, num);
}

void Profile::read(std::ifstream& fin) {
	std::string line;
	int tmp_maxl, tmp_ncodes;

	assert((fin>> tmp_maxl>> tmp_ncodes) && (tmp_ncodes == NCODES));

	if (maxL == 0) maxL = tmp_maxl;
	else assert(maxL == tmp_maxl);

	assert(p == NULL);
	allocate_space(maxL);

	for (int i = 0; i < maxL; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				assert(fin>> p[i][j][k]);

	getline(fin, line);
}

void Profile::write(std::ofstream& fout, bool isProb) {
	if (to_log_space) ss2p();

	fout<< "#pro "<< 1 + (NCODES + 1) * maxL + 1<< " Profile, format: maxL NCODES; P[POS][REF_BASE][OBSERVED_BASE], maxL blocks separated by a blank line, each block contains NCODES lines"<< std::endl;

	fout<< maxL<< '\t'<< NCODES<< std::endl;

	for (int i = 0; i < maxL; ++i) {
		for (int j = 0; j < NCODES; ++j) {
			for (int k = 0; k < NCODES - 1; ++k)
				fout<< (isProb ? p[i][j][k] : ss[i][j][k])<< '\t';
			fout<< (isProb ? p[i][j][NCODES - 1] : ss[i][j][NCODES - 1])<< std::endl;
		}
		fout<< std::endl;
	}

	fout<< std::endl;
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
				p[i][j][k] = ss[i][j][k] + pseudoC[k];
				sum += p[i][j][k];
			}
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] /= sum;
		}
}
