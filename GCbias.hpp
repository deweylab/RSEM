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

#ifndef GCBIAS_H_
#define GCBIAS_H_

#include <cassert>
#include <fstream>

#include "utils.h"

// do not estimate GCbias if nB == 0
class GCbias {
public:
	// if passed mode == INIT, do not estimate GC bias
	GCbias(model_mode_type mode, double hint = 0.04, int nB = 25);
	~GCbias();

	// return the enrichment factor for each gc bin	
	double getFactor(double gc) {
		if (mode == INIT) return 1.0;
		return gcFactor[int(gc * NUM_BIN - 1e-8)];
	}

	void udpate(double gc, double frac, bool is_forestat) {
		double *stat = (is_forestat ? forestat : backstat);
		assert(stat != NULL);
		stat[int(gc * NUM_BIN - 1e-8)] += frac;
	}

	void clear();
	void collect(const GCbias* o, bool is_forestat);
	void finish();

	void read(std::ifstream& fin, int choice); // choice: 0 -> pmf; 1 -> foreground, background
	void write(std::ofstream& fout, int choice);

private:
	static const int NUM_BIN = 100; // number of GC bins for statistics

	model_mode_type mode;
	double hint; // minimum fraction of data to become a bin
	int nB; // number of bins

	int *gcUpper; // upper bound for GC portions for each bin, [gcUpper[i - 1], gcUpper[i])

	double *forestat, *backstat, *gcFactor; // forestat, backstat: stat for each of NUM_BINs; gcFactor: factor for each NUM_BINs
	double *foreground, *background, *factor; // factor = foreground / background

	void aggregate();
	void collectSS();
};


#endif /* GCBIAS_H_ */
