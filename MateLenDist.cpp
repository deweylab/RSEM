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
#include <cassert>
#include <string>
#include <fstream>

#include "MateLenDist.hpp"

MateLenDist::MateLenDist(int mode, int maxL): mode(mode), pmf(NULL), ss(NULL), cdf(NULL) {
	lb = 1; ub = maxL; span = ub - lb + 1;
	pmf = new double[span];
	if (mode == 0) ss = new double[span];
}

MateLenDist::~MateLenDist() {
	if (pmf != NULL) delete pmf;
	if (ss != NULL) delete ss;
	if (cdf != NULL) delete cdf;
}

void findBoundaries() {
	while (ub >= lb && pmf[ub] == 0.0) --ub;
	while (lb <= ub && pmf[lb] == 0.0) ++lb;
	span = ub - lb + 1;
	assert(span > 0);
}

void MateLenDist::clear() {
	memset(pmf, 0, sizeof(double) * span);
}

void MateLenDist::collect(const MateLenDist* o) {
	for (int i = 0; i < span; ++i) pmf[i] += o->pmf[i];
}

void MateLenDist::finish() {
	memcpy(ss, pmf, sizeof(double) * span);
	ss2p();
	p2logp();
}

double MateLenDist::calcLogP() const {
	double logp = 0.0;
	for (int i = 0; i < span; ++i)
		if (ss[i] > 0.0) logp += ss[i] * pmf[i];
	return logp;
}

void MateLenDist::read(std::ifstream& fin, int choice) {
	std::string line;
	double *in = NULL;

	switch(choice) {
		case 0: in = pmf; break;
		case 1: in = ss;
	}

	assert(fin>> lb>> ub>> span);
	for (int i = 0; i < span; ++i) assert(fin>> in[i]);
	getline(fin, line);

	if (mode == 0 && choice == o) p2logp();
	if (mode == 2) prepare_for_simulation();
}

void MateLenDist::write(std::ofstream& fout, int choice) {
	double *out = NULL;

	switch(choice) {
		case 0: ss2p(); out = pmf; break;
		case 1: out = ss;
	}

	fout<< "#mld\tMateLenDist\t4\tformat: lb ub span; probability mass function values in [lb, ub], span = ub - lb + 1"<< std::endl;
	fout<< lb<< '\t'<< ub<< '\t'<< span<< std::endl;  
	for (int i = 1; i < span - 1; ++i) fout<< out[i]<< '\t';
	fout<< out[span - 1]<< std::endl<< std::endl<< std::endl;
}

void MateLenDist::ss2p() {	
	double sum;
	for (int i = 0; i < span; ++i) sum += ss[i];
	for (int i = 0; i < span; ++i) pmf[i] = ss[i] / sum;
}

void MateLenDist::p2logp() {
	for (int i = 0; i < span; ++i) pmf[i] = (pmf[i] > 0.0 : log(pmf[i]) : -std::numeric_limits<double>::infinity());
}

void MateLenDist::prepare_for_simulation() {
	if (cdf == NULL) cdf = new double[span];
	cdf[0] = pmf[0];
	for (int i = 1; i < span; ++i) cdf[i] = cdf[i - 1] + pmf[i];
}
