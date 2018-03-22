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

#ifndef FSPD_H_
#define FSPD_H_

#include <cassert>
#include <fstream>

#include "utils.h"

class FSPD {
public:
	// if mode == INIT, do not estimate FSPD
	FSPD(model_mode_type mode, int number_of_categories = 5, int number_of_bins = 10);
	~FSPD();
	
	// efflen = reflen - fragment_length + 1
	double getProb(int leftmost_pos, int efflen) {
		if (nCat == 0) return 1.0 / efflen;
		int cat = locateCategory(efflen);
		return evalCDF(leftmost_pos + 1, efflen, cat) - evalCDF(leftmost_pos, efflen, cat);
	}

	void udpate(int leftmost_pos, int efflen, double frac) {
		assert(foreground != NULL);
		int cat = locateCategory(efflen);

		int fr, to;
		double a, b;

		a = double(leftmost_pos) / efflen;
		fr = leftmost_pos * nB / efflen;
		to = ((leftmost_pos + 1) * nB - 1) / efflen;

		for (int i = fr; i < to; ++i) {
			b = double(i) / nB;
			foreground[i] += (b - a) * efflen * frac;
			a = b;
		}
		b = (leftmost_pos + 1.0) / efflen;
		foreground[i] += (b - a) * efflen * frac;
	}

	void clear();
	void collect(const FSPD* o, bool is_foreground);
	void finish();

	void read(std::ifstream& fin, int choice); // choice: 0 -> pmf; 1 -> foreground, background
	void write(std::ofstream& fout, int choice);

private:
	model_mode_type mode;
	int nCat, nB; // number of categories and number of bins

	int *catUpper; // catUpper[i] is the maximum length for category i
	double **foreground, **background, **pmf, **cdf; // pmf = foreground / background; cdf is partial sum of pmf



	int locateCategory(int efflen) {
		for (int i = 0; i < nCat; ++i) if (catUpper[i] >= efflen) return i;
		assert(false);
	}

	double evalCDF(int leftmost_pos, int efflen, int cat) {
		double val = double(leftmost_pos) * nB / efflen;
		int i = int(val - 1e-8); // assume no efflen > 1e7
	
		return cdf[cat][i] + (val - i - 1) * pmf[cat][i];
	}

	void create(double**& arr);
	void release(double** arr);
};


#endif /* FSPD_H_ */
