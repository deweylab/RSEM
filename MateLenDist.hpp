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

#ifndef MATELENDIST_H_
#define MATELENDIST_H_

#include <cassert>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class MateLenDist {
public:
	MateLenDist(model_mode_type mode, int lb, int ub);
	~MateLenDist();

	void findBoundaries();

	int getMinL() const { return lb; }
	int getMaxL() const { return ub; }
	
	double getLogProb(int len) const { 
		assert(len <= ub);
		return (len >= lb ? pmf[len - lb] : pmf[0]);
	}
	
	void update(int len) { 
		assert(len <= ub);
		if (len >= lb) ++pmf[len - lb];
	}
	
	void clear();
	void collect(const MateLenDist* o);
	void finish();
	
	double calcLogP() const;
	
	void read(std::ifstream& fin, int choice); // choice: 0 -> p; 1 -> ss
	void write(std::ofstream& fout, int choice);
	
	int simulate(Sampler* sampler) {
		return lb + sampler->sample(cdf, span);
	}

private:
	model_mode_type mode;
	int lb, ub, span; // [lb, ub], span = ub - lb + 1
	double *pmf, *ss, *cdf; // probability mass function, sufficiant statistics, and cumulative density function

	void prepare_for_simulation();

	void ss2p(); // from sufficient statistics to pmf
	void p2logp(); // convert to log space
};

#endif /* MATELENDIST_H_ */
