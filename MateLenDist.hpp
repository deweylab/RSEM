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

#ifndef MATELENDIST_H_
#define MATELENDIST_H_

#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>

#include "sampling.hpp"

class MateLenDist {
public:
	MateLenDist();
	
	int getMinL() const { return lb; }
	int getMaxL() const { return ub; }

	/*
		comment: In theory, this length should be the untrimmed length. If necessary, may ask user to provide the untrimmed lengths of each trimmed read (for example, 
						 encode this information in the header line of each read).
						 Do not need this for the EM algorithm, only need it for calculating the complete log likelihood score.
	*/
	double getLogProb(int len) const {
		assert(len > lb && len <= ub);
		return log(pmf[len - lb]);
	}
	
	void update(int len) {
		if (len > ub) { ub = len; ss.resize(ub, 0.0); }
		++ss[len - lb];
	}
	
	void finish();
	
	/*
		comment: Calculate log probability for unalignable reads, call after finish()
	*/
	double calcLogP() const;
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, int mate, bool isProb = true); // mate = 1 or 2

	void prepare_for_simulation();
	
	// In the RNASeqModel, if the simulate mate length > fragment length, the fragment length is used
	int simulate(Sampler* sampler) {
		return lb + sampler->sample(cdf, span);
	}

private:
	int lb, ub, span; // [lb, ub], span = ub - lb + 1
	std::vector<double> pmf, cdf, ss; // probability mass function, cumulative density function and sufficiant statistics
};

#endif 
