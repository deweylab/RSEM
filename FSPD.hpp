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
	// if passed mode == INIT, no estimation of FSPD
	FSPD(model_mode_type mode, double hint = 0.2, int nCat = 5, int nB = 10);
	~FSPD();
	
	// efflen = reflen - fragment_length + 1
	double getProb(int leftmost_pos, int efflen) {
		if (mode == INIT) return 1.0 / efflen;
		int cat = catMap[(efflen < MAXL_FSPD ? efflen / INTERVAL : NUM_CAT - 1)];
		return evalCDF(leftmost_pos + 1, efflen, cat) - evalCDF(leftmost_pos, efflen, cat);
	}

	void udpate(int leftmost_pos, int efflen, double frac, bool is_forestat) {
		int cat = (efflen < MAXL_FSPD ? efflen / INTERVAL : NUM_CAT - 1);
		double **stat = (is_forestat ? forestat : backstat);
		int fr, to;
		double a, b;

		assert(stat != NULL);

		a = double(leftmost_pos) / efflen;
		fr = leftmost_pos * nB / efflen;
		to = ((leftmost_pos + 1) * nB - 1) / efflen;

		for (int i = fr; i < to; ++i) {
			b = double(i) / nB;
			stat[cat][i] += (b - a) * efflen * frac;
			a = b;
		}
		b = (leftmost_pos + 1.0) / efflen;
		stat[cat][i] += (b - a) * efflen * frac;
	}

	void clear();
	void collect(const FSPD* o, bool is_forestat);
	void finish();

	void read(std::ifstream& fin, int choice); // choice: 0 -> pmf; 1 -> foreground, background
	void write(std::ofstream& fout, int choice);

private:
	static const int MAXL_FSPD = 10000; // maximum length in FSPD category, if l >= 10000, add the count to bin of 10000
	static const int INTERVAL = 100; // store statistics every 100bp
	static const int NUM_CAT = 100; // number of categories for statistics

	model_mode_type mode;
	double hint; // minimum fraction of data to become a separate category
	int nCat, nB; // number of categories and number of bins

	int *catUpper; // In NUM_CAT intervals, category i contains intervals [catUpper[i - 1], catUpper[i]) 
	int *catMap; // category map, from NUM_CAT to naCat

	double **forestat, **backstat; // forestat, backstat, NUM_CAT x nB, [i * INTERVAL, (i + 1) * INTERVAL)
	double **foreground, **background, **pmf, **cdf; // pmf = foreground / background; cdf is partial sum of pmf

	double evalCDF(int leftmost_pos, int efflen, int cat) {
		double val = double(leftmost_pos) * nB / efflen;
		int i = int(val - 1e-8); // assume no efflen > 1e7
	
		return cdf[cat][i] + (val - i - 1) * pmf[cat][i];
	}

	void create(double**& arr, int num);
	void release(double** arr, int num);
	void aggregate(); // from refined NUM_CAT to nCat, automatically adjust cat size based on data
};


#endif /* FSPD_H_ */
