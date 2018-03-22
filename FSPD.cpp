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
#include "FSPD.hpp"

FSPD::FSPD(model_mode_type mode, int number_of_categories, int number_of_bins) : mode(mode), nCat(number_of_categories), nB(number_of_bins) {
	catUpper = NULL;
	foreground = background = pmf = cdf = NULL;

	if (mode == INIT) { nCat = nB = 0; }
	if (nCat == 0) return;
	if (nCat > 1) catUpper = new int[nCat - 1];
	create(pmf);
	create(cdf);
	if (mode == MASTER || mode == CHILD) {
		create(foreground);
		create(background);
	}
}

~FSPD::FSPD() {
	if (catUpper != NULL) delete[] catUpper;
	if (foreground != NULL) release(foreground);
	if (background != NULL) release(background);
	if (pmf != NULL) release(pmf);
	if (cdf != NULL) release(cdf);
}

void FSPD::calcCatUpper(const double* weights) {
	int pos;
	double sum = 0.0, one_slice, expected;

	for (int i = 1; i <= MAXL_FSPD; ++i) sum += weights[i];
	one_slice = sum / nCat;

	pos = 0;
	expected = sum = 0.0;
	for (int i = 0; i < nCat - 1; ++i) {
		expected += one_slice;
		do {
			sum += weights[++pos]; 
		} while (sum < expected && pos <= MAXL_FSPD - (nCat - 1 - i));
		if (sum >= expected && (expected - (sum - weights[pos])) < (sum - expected)) sum -= weights[pos--];
		catUpper[i] = pos;
	}
}

void FSPD::clear() {
	for (int i = 0; i < nCat; ++i) {
		memset(foreground[i], 0, sizeof(double) * nB);
		memset(background[i], 0, sizeof(double) * nB);
	}
}

void FSPD::collect(const FSPD* o, bool is_foreground) {
	for (int i = 0; i < nCat; ++i)
		for (int j = 0; j < nB; ++j)
			if (is_foreground) foreground[i][j] += o->foreground[i][j];
			else background[i][j] += o->background[i][j];
}

void FSPD::finish() {
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
	int tmp_ncat, tmp_nb;

	assert((fin>> line) && (line == "#fspd"));
	assert(getline(fin, line));
	assert((fin>> tmp_ncat>> tmp_nb) && (tmp_ncat == nCat) && (tmp_nb == nB));
	if (nCat == 0) return;

	for (int i = 0; i < nCat; ++i) assert(fin>> catUpper[i]);		

	if (choice == 0) {
		for (int i = 0; i < nCat; ++i) 
			for (int j = 0; j < nB; ++i) {
				assert(fin>> pmf[i][j]);
				cdf[i][j] = pmf[i][j];
				if (j > 0) cdf[i][j] += cdf[i][j - 1];
			}
	}
	else {
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
	if (nCat == 0) {
		fout<< "#fspd\tFragment Start Position Distribution\t3\tformat: nCat nB, nCat == 0 means no FSPD estimation"<< std::endl;
		fout<< nCat<< "\t"<< nB<< std::endl<< std::endl<< std::endl;
		return;
	}

	switch(choice) {
		case 0:
			fout<< "#fspd\tFragment Start Position Distribution\t"<< nCat + 5<< "\tformat: nCat nB; nCat values, category upper bound; nCat x nB values for pmf"<< std::endl;
			fout<< nCat<< '\t'<< nB<< std::endl;
			for (int i = 0; i < nCat - 1; ++i) fout<< catUpper[i]<< '\t';
			fout<< catUpper[nCat - 1]<< std::endl<< std::endl;
			for (int i = 0; i < nCat; ++i) {
				for (int j = 0; j < nB - 1; ++j) fout<< pmf[i][j]<< '\t';
				fout<< pmf[i][nB - 1]<< std::endl;
			}
			fout<< std::endl<< std::endl;
			break;
		case 1:
			fout<< "#fspd\tFragment Start Position Distribution\t"<< nCat * 2 + 6<< "\tformat: nCat nB; nCat values, category upper bound; nCat x nB values for foreground; nCat x nB values for background"<< std::endl;
			fout<< nCat<< '\t'<< nB<< std::endl;
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

void FSPD::create(double**& arr) {
	arr = new double*[nCat];
	for (int i = 0; i < nCat; ++i) arr[i] = new double[nB];
}

void FSPD::release(double** arr) {
	for (int i = 0; i < nCat; ++i) delete[] arr[i];
	delete[] arr;
}
