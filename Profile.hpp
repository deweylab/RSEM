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

#ifndef PROFILE_H_
#define PROFILE_H_

#include <cstdlib>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Profile {
public:
	Profile(bool to_log_space = false, int maxL = 0);
	~Profile();

	double getLogProb(int pos, int ref_base, int read_base) const {
		return p[pos][ref_base][read_base]; // p is in log space
	}

	void expand(int len) {
		if (len <= maxL) return;
		p = (double (**)[NCODES])realloc(p, sizeof(double (*)[NCODES]) * len);
		for (int i = maxL; i < len; ++i)
			p[i] = (double (*)[NCODES])calloc(NCODES * NCODES, sizeof(double));
		maxL = len;
	}

	void update(int pos, int ref_base, int read_base, double frac) {
		p[pos][ref_base][read_base] += frac;
	}

	void collect(const Profile* o);
	void finish(int length = -1);
	void clear();
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, bool isProb = true);

	void prepare_for_simulation();
	
	char simulate(Sampler* sampler, int pos, int ref_base) {
		return code2base[sampler->sample(p[pos][ref_base], NCODES)];
	}
	
private:
	bool to_log_space; // true if we should transfer p to log space

	int maxL; // profile length
	double (**p)[NCODES]; //profile matrices
	double (*ss)[NCODES][NCODES];

	void ss2p(); // calculate p from ss

	void allocate_space(int len) {
		p = (double (**)[NCODES])calloc(len, sizeof(double (*)[NCODES]));
		for (int i = 0; i < len; ++i)
			p[i] = (double (*)[NCODES])calloc(NCODES * NCODES, sizeof(double)); // should be 0 initialized
	}
};

#endif /* PROFILE_H_ */
