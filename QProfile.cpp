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
#include "QProfile.hpp"

QProfile::QProfile(bool to_log_space) : to_log_space(to_log_space) {
	memset(p, 0, sizeof(p));	
	ss = NULL;
}

QProfile::~QProfile() {
	if (ss != NULL) delete[] ss;
}

void QProfile::collect(const QProfile* o) {
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				p[i][j][k] += o->p[i][j][k];
}

void QProfile::finish() {
	// collect sufficient statistics
	if (ss == NULL) ss = new double[QSIZE][NCODES][NCODES];
	memcpy(ss, p, sizeof(p));
	// calculate p
	ss2p();
	// convert to log space
	if (to_log_space)
		for (int i = 0; i < QSIZE; ++i)
			for (int j = 0; j < NCODES; ++j)
				for (int k = 0; k < NCODES; ++k)
					p[i][j][k] = log(p[i][j][k]);
}

void QProfile::clear() {
	memset(p, 0, sizeof(p));
}

void QProfile::read(std::ifstream& fin) {
	std::string line;
	int tmp_qsize, tmp_ncodes;

	assert((fin>> tmp_qsize>> tmp_ncodes) && (tmp_qsize == QSIZE) && (tmp_ncodes == NCODES));
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 0; k < NCODES; ++k)
				assert(fin>> p[i][j][k]);
	getline(fin, line);
}

void QProfile::write(std::ofstream& fout, bool isProb) {
	if (isProb && to_log_space) ss2p();

	fout<< "#qpro "<< 1 + (NCODES + 1) * QSIZE + 1<< " QProfile, format: SIZE NCODES; P[QUAL][REF_BASE][OBSERVED_BASE], SIZE blocks separated by a blank line, each block contains NCODES lines"<< std::endl;

	fout<< QSIZE<< '\t'<< NCODES<< std::endl;
	for (int i = 0; i < QSIZE; ++i) {
		for (int j = 0; j < NCODES; ++j) {
			for (int k = 0; k < NCODES - 1; ++k)
				fout<< (isProb ? p[i][j][k] : ss[i][j][k])<< '\t';
			fout<< (isProb ? p[i][j][NCODES - 1] : ss[i][j][NCODES - 1])<< std::endl;
		}
		fout<< std::endl;
	}
	fout<< std::endl;
}

void QProfile::prepare_for_simulation() {
	for (int i = 0; i < QSIZE; ++i)
		for (int j = 0; j < NCODES; ++j)
			for (int k = 1; k < NCODES; ++k)
				p[i][j][k] += p[i][j][k - 1];
}

void QProfile::ss2p() {
	double sum;
	
	for (int i = 0; i < QSIZE; ++i) 
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
