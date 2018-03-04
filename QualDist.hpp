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

#ifndef QUALDIST_H_
#define QUALDIST_H_

#include <cmath>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"
#include "QUALstring.hpp"

//from 33 to 126 to encode 0 to 93
class QualDist {
public:
	QualDist();
	~QualDist();
	
	void update(const QUALstring& qual) {
		int len = qual.getLen();

		++p_init[qual.qualAt(0)];
		for (int i = 1; i < len; ++i) {
			++p_tran[qual.qualAt(i - 1)][qual.qualAt(i)];
		}
	}

	void collect(const QualDist* o);
	void finish();
	
	double getLogProb(const QUALstring& qual) const {
		int len = qual.getLen();
		double log_prob = 0.0;
		
		log_prob += log(p_init[qual.qualAt(0)]);
		for (int i = 1; i < len; ++i) {
			log_prob += log(p_tran[qual.qualAt(i - 1)][qual.qualAt(i)]);
		}
		
		return log_prob;
	}

	double calcLogP() const;
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, bool isProb = true);
	
	void prepare_for_simulation();

	// The first two simulate functions are used for simulating alignable reads
	int simulate(Sampler* sampler) {
		return sampler->sample(p_init, QSIZE);
	}

	int simulate(Sampler* sampler, int old_qval) {
		return sampler->sample(p_tran[old_qval], QSIZE);
	}

	void simulate(Sampler* sampler, int len, std::string& qual) {
		int qval;

		qual.assign(len, 0);
		qval = sampler->sample(p_init, QSIZE);
		qual[0] = qval2char(qval);
		for (int i = 1; i < len; ++i) {
			qval = sampler->sample(p_tran[qval], QSIZE);
			qual[i] = qval2char(qval);
		}
	}
	
private:
	double p_init[QSIZE];
	double p_tran[QSIZE][QSIZE]; //p_tran[a][b] = p(b|a)

	double *ss_init;
	double (*ss_tran)[QSIZE];
};

#endif /* QUALDIST_H_ */
