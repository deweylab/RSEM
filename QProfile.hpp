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

#ifndef QPROFILE_H_
#define QPROFILE_H_

#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class QProfile {
public:
	QProfile(model_mode_type mode);
	~QProfile();
	
	// qual starts from 0, 33 is already deducted
	double getLogProb(int qual, int ref_base, int read_base) const {
		return p[qual][ref_base][read_base]; // p in log space
	}

	void update(int qual, int ref_base, int read_base, double frac) {
		p[qual][ref_base][read_base] += frac;
	}

	void clear();
	void collect(const QProfile* o);
	void finish();
	
	void read(std::ifstream& fin, int choice); // choice: 0 -> p; 1 -> ss; 2 -> c; 
	void write(std::ofstream& fout, int choice);

	char simulate(Sampler* sampler, int qual, int ref_base) {
		return code2base[sampler->sample(p[qual][ref_base], NCODES)];
	}
	
private:
	model_mode_type mode;
	double (*p)[NCODES][NCODES], (*ss)[NCODES][NCODES]; // p[q][r][c] = p(c|r,q), if mode == 0, p is in log space; ss, sufficient statistics

	void init();
	void prepare_for_simulation();

	void ss2p(); // from sufficient statistics to p
	void p2logp(); // convert to log space
};

#endif /* QPROFILE_H_ */
