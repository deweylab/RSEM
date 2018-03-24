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
#include <algorithm>

#include "utils.h"
#include "FSPD.hpp"

FSPD::FSPD(model_mode_type mode, double hint, int nCat, int nB) : mode(mode), hint(hint), nCat(nCat), nB(nB) {
 	catUpper = catMap = NULL;
	forestat = backstat = NULL;
	foreground = background = pmf = cdf = NULL;

	if (mode == INIT) return;

	int maxCat = (mode == SIMULATION ? nCat : int(1.0 / hint + 1e-8));

	if (mode != CHILD) { catUpper = new int[maxCat]; catMap = new int[NUM_CAT]; create(pmf, maxCat); create(cdf, maxCat); }
	if (mode == FIRST_PASS || mode == MASTER) { create(foreground, maxCat); create(background, maxCat); }
	if (mode == MASTER || mode == CHILD) { create(forestat, NUM_CAT); create(backstat, NUM_CAT); }
}

~FSPD::FSPD() {
	if (catUpper != NULL) delete[] catUpper;
	if (catMap != NULL) delete[] catMap;
	if (forestat != NULL) release(forestat, NUM_CAT);
	if (backstat != NULL) release(backstat, NUM_CAT);
	if (foreground != NULL) release(foreground, nCat);
	if (background != NULL) release(background, nCat);
	if (pmf != NULL) release(pmf, nCat);
	if (cdf != NULL) release(cdf, nCat);
}

void FSPD::clear() {
	for (int i = 0; i < NUM_CAT; ++i) {
		memset(forestat[i], 0, sizeof(double) * nB);
		memset(backstat[i], 0, sizeof(double) * nB);
	}
}

void FSPD::collect(const FSPD* o, bool is_forestat) {
	for (int i = 0; i < NUM_CAT; ++i)
		for (int j = 0; j < nB; ++j)
			if (is_forestat) forestat[i][j] += o->forestat[i][j];
			else backstat[i][j] += o->backstat[i][j];
}

void FSPD::finish() {
	if (mode == FIRST_PASS || mode == MASTER) aggregate(); // could change it to only FIRST_PASS

	// compute catMap
	int fr = 0;
	for (int i = 0; i < nCat; ++i) 
		while (fr < catUpper[i]) catMap[fr++] = i;
	// calculate pmf and cdf
	for (int i = 0; i < nCat; ++i) {
		for (int j = 0; j < nB; ++j) {
			pmf[i][j] = (foreground[i][j] + pseudo_count_FSPD) / (background[i][j] + pseudo_count_FSPD);
			cdf[i][j] = pmf[i][j];
			if (j > 0) cdf[i][j] += cdf[i][j - 1];
		}
		// normalization
		for (int j = 0; j < nB; ++j) {
			pmf[i][j] /= cdf[i][nB - 1];
			cdf[i][j] /= cdf[i][nB - 1];
		}
	}
}

void FSPD::read(std::ifstream& fin, int choice) {
	std::string line;
	double tmp_hint;
	int tmp_ncat, tmp_nb;

	assert((fin>> line) && (line == "#fspd"));
	assert(getline(fin, line));
	if (mode == INIT) return;

	assert((fin>> tmp_ncat>> tmp_nb>> tmp_hint) && (tmp_ncat == nCat) && (tmp_nb == nB) && (tmp_hint == hint));
	for (int i = 0; i < nCat; ++i) assert(fin>> catUpper[i]);

	switch(choice) {
		case 0:
			int fr =0 ;
			for (int i = 0; i < nCat; ++i) {
				while (fr < catUpper[i]) catMap[fr++] = i;
				
				for (int j = 0; j < nB; ++j) {
					assert(fin>> pmf[i][j]);
					cdf[i][j] = pmf[i][j];
					if (j > 0) cdf[i][j] += cdf[i][j - 1];
				}
			}
			break;
		case 1:
			for (int i = 0; i < nCat; ++i)
				for (int j = 0; j < nB; ++j)
					assert(fin>> foreground[i][j]);
			for (int i = 0; i < nCat; ++i)
				for (int j = 0; j < nB; ++j)
					assert(fin>> background[i][j]);
	}
	assert(getline(fin, line));
}

void FSPD::write(std::ofstream& fout, int choice) {
	if (mode == INIT) {
		fout<< "#fspd\tFragment Start Position Distribution\t2\tno FSPD estimation"<< std::endl<< std::endl<< std::endl;
		return;
	}

	switch(choice) {
		case 0:
			fout<< "#fspd\tFragment Start Position Distribution\t"<< nCat + 5<< "\tformat: nCat nB hint; nCat values, category upper bound; nCat x nB values for pmf"<< std::endl;
			fout<< nCat<< '\t'<< nB<< '\t'<< hint<< std::endl;
			for (int i = 0; i < nCat - 1; ++i) fout<< catUpper[i]<< '\t';
			fout<< catUpper[nCat - 1]<< std::endl<< std::endl;
			for (int i = 0; i < nCat; ++i) {
				for (int j = 0; j < nB - 1; ++j) fout<< pmf[i][j]<< '\t';
				fout<< pmf[i][nB - 1]<< std::endl;
			}
			fout<< std::endl<< std::endl;
			break;
		case 1:
			fout<< "#fspd\tFragment Start Position Distribution\t"<< nCat * 2 + 6<< "\tformat: nCat nB hint; nCat values, category upper bound; nCat x nB values for foreground; nCat x nB values for background"<< std::endl;
			fout<< nCat<< '\t'<< nB<< '\t'<< hint<< std::endl;
			for (int i = 0; i < nCat - 1; ++i) fout<< catUpper[i]<< '\t';
			fout<< catUpper[nCat - 1]<< std::endl<< std::endl;
			for (int i = 0; i < nCat; ++i) {
				for (int j = 0; j < nB - 1; ++j) fout<< foreground[i][j]<< '\t';
				fout<< foreground[i][nB - 1]<< std::endl;
			}
			fout<< std::endl;
			for (int i = 0; i < nCat; ++i) {
				for (int j = 0; j < nB - 1; ++j) fout<< background[i][j]<< '\t';
				fout<< background[i][nB - 1]<< std::endl;
			}
			fout<< std::endl<< std::endl;
	}
}

void FSPD::create(double**& arr, int num) {
	arr = new double*[num];
	for (int i = 0; i < num; ++i) arr[i] = new double[nB];
}

void FSPD::release(double** arr, int num) {
	for (int i = 0; i < num; ++i) delete[] arr[i];
	delete[] arr;
}

void FSPD::aggregate() {
	int pos, nzero;
	double psum[NUM_CAT];
	double sum, one_slice;

	sum = 0.0;
	for (int i = 0; i < NUM_CAT; ++i) {
		psum[i] = 0.0;
		for (int j = 0; j < nB; ++j) psum[i] += backstat[i][j];
		sum += psum[i];
	}

	one_slice = sum * hint;	
	nCat = 0; pos = 0;
	while (pos < NUM_CAT) {
		sum = 0.0;
		do { 
			sum += psum[pos++];
		} while (sum < one_slice && pos < NUM_CAT);
		nzero = 0;
		while (pos < NUM_CAT && psum[pos] == 0.0) ++pos, ++nzero;
		catUpper[nCat++] = pos - nzero / 2;
	}
	// if the last category's mass is less than 0.8 one_slice, make it independent, otherwise merge it into the one on the left
	if (sum < one_slice * 0.8) catUpper[(--nCat) - 1] = pos;
}

// fill in foreground and background
void FSPD::collectSS() {
	int pos = 0;
	for (int i = 0; i < nCat; ++i) {
		memset(foreground[i], 0, sizeof(double) * nB);
		memset(background[i], 0, sizeof(double) * nB);
		while (pos < catUpper[i]) {
			for (int j = 0; j < nB; ++j) {
				foreground[i][j] += forestat[pos][j];
				background[i][j] += backstat[pos][j];
			}
			++pos;
		}
	}
}
