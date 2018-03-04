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

#ifndef QPROFILE_H_
#define QPROFILE_H_

#include <cassert>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class QProfile {
public:
	QProfile(bool to_log_space = false);
	~QProfile();
	
	// qual starts from 0, 33 is already deducted
	double getLogProb(int qual, int ref_base, int read_base) const {
		return p[qual][ref_base][read_base]; // p in log space
	}

	void update(int qual, int ref_base, int read_base, double frac) {
		p[qual][ref_base][read_base] += frac;
	}

	void collect(const QProfile* o);
	void finish();
	void clear();
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, bool isProb);

	void prepare_for_simulation();
	
	char simulate(Sampler* sampler, int qual, int ref_base) {
		return code2base[sampler->sample(p[qual][ref_base], NCODES)];
	}
	
private:
	bool to_log_space; // if we should transfer p to log space
	
	double p[QSIZE][NCODES][NCODES]; // p[q][r][c] = p(c|r,q), if isMaster == true, p is in log space
	double (*ss)[NCODES][NCODES]; // sufficient statistics

	void ss2p(); // from sufficient statistics to p
};

#endif /* QPROFILE_H_ */
