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
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "GCbias.hpp"

GCbias::GCbias(model_mode_type mode, double hint, int nB) : mode(mode), hint(hint), nB(nB) {
 	gcUpper = NULL;
	forestat = backstat = gcFactor = NULL;
	foreground = background = factor = NULL;

	if (mode == INIT) return;

	int maxB = (mode == SIMULATION ? nB : int(1.0 / hint + 1e-8));

	if (mode != CHILD) { gcUpper = new int[maxB]; gcMap = new int[NUM_BIN]; factor = new double[maxB]; }
	if (mode == FIRST_PASS || mode == MASTER) { foreground = new double[maxB]; background = new double[maxB]; }
	if (mode == MASTER || mode == CHILD) { forestat = new double[NUM_BIN]; backstat = new double[NUM_BIN]; }
}

~GCbias::GCbias() {
	if (gcUpper != NULL) delete[] gcUpper;
	if (forestat != NULL) delete[] forestat;
	if (backstat != NULL) delete[] backstat;
	if (gcFactor != NULL) delete[] gcFactor;
	if (foreground != NULL) delete[] foregorund;
	if (background != NULL) delete[] background;
	if (factor != NULL) delete[] factor;
}

void GCbias::clear() {
	memset(forestat, 0, sizeof(double) * NUM_BIN);
	memset(backstat, 0, sizeof(double) * NUM_BIN);
}

void GCbias::collect(const GCbias* o, bool is_forestat) {
	for (int i = 0; i < NUM_BIN; ++i)
		if (is_forestat) forestat[i] += o->forestat[i];
		else backstat[i] += o->backstat[i];
}

void GCbias::finish() {
	double foresum, backsum, denom;

	if (mode == FIRST_PASS || mode == MASTER) aggregate();
	collectSS();

	foresum = backsum = 0.0;
	for (int i = 0; i < nB; ++i) {
		foresum += foreground[i] + pseudo_count_GCbias;
		backsum += background[i] + pseudo_count_GCbias;
		factor[i] = (foreground[i] + pseudo_count_GCbias) / (background[i] + pseudo_count_GCbias);
	}
	denom = foresum / backsum;
	for (int i = 0; i < nB; ++i) factor[i] /= denom;

	int fr = 0;
	for (int i = 0; i < nB: ++i)
		while (fr < gcUpper[i]) gcFactor[fr++] = factor[i];
}

void GCbias::read(std::ifstream& fin, int choice) {
	std::string line;
	double tmp_hint;
	int tmp_nb;

	assert((fin>> line) && (line == "#gcbias"));
	assert(getline(fin, line));
	if (mode == INIT) return;

	assert((fin>> tmp_nb>> tmp_hint) && (tmp_nb == nB) && (tmp_hint == hint));
	for (int i = 0; i < nB; ++i) assert(fin>> gcUpper[i]);

	switch(choice) {
		case 0: 
			int fr = 0;
			for (int i = 0; i < nB; ++i) {
				assert(fin>> factor[i]);
				while (fr < gcUpper[i]) gcFactor[fr++] = factor[i];
			}
			break;
		case 1:
			for (int i = 0; i < nB; ++i) assert(fin>> foreground[i]);
			for (int i = 0; i < nB; ++i) assert(fin>> background[i]);
	}
	assert(getline(fin, line));
}

void GCbias::write(std::ofstream& fout, int choice) {
	if (mode == INIT) {
		fout<< "#gcbias\tGC content bias model\t2\tformat: no GC bias estimation"<< std::endl<< std::endl<< std::endl;
		return;
	}

	switch(choice) {
		case 0:
			fout<< "#gcbias\tGC content bias model\t5\tformat: nB hint; nB values, gc upper bound; nB values, factors"<< std::endl;
			fout<< nB<< '\t'<< hint<< std::endl;
			for (int i = 0; i < nB - 1; ++i) fout<< gcUpper[i]<< '\t';
			fout<< gcUpper[nB - 1]<< std::endl;
			for (int i = 0; i < nB - 1; ++i) fout<< factor[i]<< '\t';
			fout<< factor[nB - 1]<< std::endl<< std::endl<< std::endl;
			break;
		case 1:
			fout<< "#gcbias\tGC content bias model\t6\tformat: nB hint; nB values, gc upper bound; nB values, foreground; nB values, background"<< std::endl;
			fout<< nB<< '\t'<< hint<< std::endl;
			for (int i = 0; i < nB - 1; ++i) fout<< gcUpper[i]<< '\t';
			fout<< gcUpper[nB - 1]<< std::endl;
			for (int i = 0; i < nB - 1; ++i) fout<< foreground[i]<< '\t';
			fout<< foreground[nB - 1]<< std::endl;
			for (int i = 0; i < nB - 1; ++i) fout<< background[i]<< '\t';
			fout<< background[nB - 1]<< std::endl<< std::endl<< std::endl;
	}
}

void GCbias::aggregate() {
	int pos, nzero;
	double sum, one_slice;

	sum = 0.0;
	for (int i = 0; i < NUM_BIN; ++i) sum += backstat[i];

	one_slice = sum * hint;
	nB = 0; pos = 0;
	while (pos < NUM_BIN) {
		sum = 0.0;
		do {
			sum += bacstat[pos++];
		} while (sum < one_slice && pos < NUM_BIN);
		nzero =0 ;
		while (pos < NUM_BIN && backstat[pos] == 0.0) ++pos, ++nzero;
		gcUpper[nB++] = pos - nzero / 2;
	}

	if (sum < one_slice * 0.8) gcUpper[(--nB) - 1] = pos;
}

// collect sufficient statistics
void GCbias::collectSS() {
	int	pos = 0;
	for (int i = 0; i < nB; ++i) {
		foreground[i] = background[i] = 0.0;
		while (pos < gcUpper[i]) {
			foreground[i] += forestat[pos];
			background[i] += backstat[pos];
			++pos;
		}
	}
}
