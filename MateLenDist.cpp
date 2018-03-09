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

#include "utils.h"
#include "MateLenDist.hpp"

MateLenDist::MateLenDist(model_mode_type mode, int lb, int ub): mode(mode), lb(lb), ub(ub), span(ub - lb + 1), pmf(NULL), ss(NULL), cdf(NULL) {
	pmf = new double[span];
	if (mode == FIRST_PASS || mode == MASTER) ss = new double[span];
	if (mode == SIMULATION) cdf = new double[span];
}

MateLenDist::~MateLenDist() {
	if (pmf != NULL) delete pmf;
	if (ss != NULL) delete ss;
	if (cdf != NULL) delete cdf;
}

void MateLenDist::findBoundaries() {
	while (span > 0 && pmf[span - 1] == 0.0) --span;
	assert(span > 0);
	ub = lb + span - 1;

	int pos = 0;
	while (pos < span && pmf[pos] == 0.0) ++pos;
	lb += pos;
	span -= pos;

	if (pos > 0) memmove(pmf, pmf + pos, sizeof(double) * span);
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
	if (mode == MASTER) p2logp();
}

double MateLenDist::calcLogP() const {
	double logp = 0.0;
	for (int i = 0; i < span; ++i)
		if (ss[i] > 0.0) logp += ss[i] * pmf[i];
	return logp;
}

void MateLenDist::read(std::ifstream& fin, int choice) {
	std::string line;
	int tmp_lb, tmp_ub, tmp_span;
	double *in = NULL;

	switch(choice) {
		case 0: in = pmf; break;
		case 1: in = ss;
	}

	assert((fin>> line) && (line == "#mld"));
	assert(getline(fin, line));
	assert((fin>> tmp_lb>> tmp_ub>> tmp_span) && (tmp_lb == lb) && (tmp_ub == ub) && (tmp_span == span));
	for (int i = 0; i < span; ++i) assert(fin>> in[i]);
	assert(getline(fin, line));

	if (mode == MASTER && choice == 0) p2logp();
	if (mode == SIMULATION) prepare_for_simulation();
}

void MateLenDist::write(std::ofstream& fout, int choice) {
	double *out = NULL;

	switch(choice) {
		case 0: if (mode == MASTER) ss2p(); out = pmf; break;
		case 1: out = ss;
	}

	fout<< "#mld\tMateLenDist\t4\tformat: lb ub span; probability mass function values in [lb, ub], span = ub - lb + 1"<< std::endl;
	fout<< lb<< '\t'<< ub<< '\t'<< span<< std::endl;  
	for (int i = 1; i < span - 1; ++i) fout<< out[i]<< '\t';
	fout<< out[span - 1]<< std::endl<< std::endl<< std::endl;
}

void MateLenDist::prepare_for_simulation() {
	cdf[0] = pmf[0];
	for (int i = 1; i < span; ++i) cdf[i] = cdf[i - 1] + pmf[i];
}

void MateLenDist::ss2p() {	
	double sum = 0.0;
	for (int i = 0; i < span; ++i) sum += ss[i];
	for (int i = 0; i < span; ++i) pmf[i] = ss[i] / sum;
}

void MateLenDist::p2logp() {
	for (int i = 0; i < span; ++i) pmf[i] = (pmf[i] > 0.0 ? log(pmf[i]) : -std::numeric_limits<double>::infinity());
}
