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

#ifndef NOISEQPROFILE_H_
#define NOISEQPROFILE_H_

#include <string>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"

class NoiseQProfile {
public:
	NoiseQProfile(int mode); 
	~NoiseQProfile();

	void update(const QUALstring& qual, const SEQstring& readseq, bool is_aligned, double frac = 1.0) {
		int len = readseq.getLen();
		for (int i = 0; i < len; ++i)
			if (is_aligned) p[qual.qualAt(i)][readseq.baseCodeAt(i)] += frac;
			else c[qual.qualAt(i)][readseq.baseCodeAt(i)] += frac;
	}

	void clear();
	void collect(const NoiseQProfile* o);
	void finish();

	double getLogProb(const QUALstring& qual, const SEQstring& readseq) const {
		int len = readseq.getLen();
		double log_prob = 0.0;

		for (int i = 0; i < len; ++i)
			log_prob += p[qual.qualAt(i)][readseq.baseCodeAt(i)]; // p in log space

		return log_prob;
	}

	double calcLogP() const; // for noise reads

	void read(std::ifstream& fin, int choice); // choice: 0 -> p; 1 -> ss; 2 -> c; 
	void write(std::ofstream& fout, int choice);

	// here p is not in log space
	void simulate(Sampler *sampler, int len, const std::string& qual, std::string& readseq) {
		readseq.assign(len, 0);
		for (int i = 0; i < len; ++i)
			readseq[i] = code2base[sampler->sample(p[char2qval(qual[i])], NCODES)];
	}

private:
	int mode; // 0, master; 1, child; 2, simulation
	double (*p)[NCODES], (*c)[NCODES], (*ss)[NCODES]; // p, probability in log space; c, counts from unaligned reads; ss, sufficient statistics	

	void calc_ss(); // calculate sufficient statistics
	void ss2p(); // from sufficient statistics to probabilities
	void p2logp(); // convert to log space
	void prepare_for_simulation();
};

#endif /* NOISEQPROFILE_H_ */
