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

#ifndef PROFILE_H_
#define PROFILE_H_

#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Profile {
public:
	Profile(int mode, int maxL);
	~Profile();

	double getLogProb(int pos, int ref_base, int read_base) const {
		return p[pos][ref_base][read_base]; // p is in log space
	}

	void update(int pos, int ref_base, int read_base, double frac) {
		p[pos][ref_base][read_base] += frac;
	}

	void setMaxL(int maxL) { this->maxL = maxL; } // set maxL, only used in parseAlignments

	void clear();
	void collect(const Profile* o);
	void finish();
	
	void read(std::ifstream& fin, int choice); // choice: 0 -> p; 1 -> ss
	void write(std::ofstream& fout, int choice); 
	
	char simulate(Sampler* sampler, int pos, int ref_base) {
		return code2base[sampler->sample(p[pos][ref_base], NCODES)];
	}
	
private:
	int mode; // 0, master; 1, child; 2, simulation
	int maxL; // profile length
	double (*p)[NCODES][NCODES], (*ss)[NCODES][NCODES]; // p, probability in log space; ss, sufficient statistics

	void ss2p(); // from sufficient statistics to probabilities
	void p2logp(); // convert to log space
	void prepare_for_simulation();
};

#endif /* PROFILE_H_ */
