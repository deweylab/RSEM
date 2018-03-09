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

#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Orientation {
 public:
	Orientation(model_mode_type mode, double probF = 0.5);
	~Orientation();

	//dir : +/-
	double getProb(char dir) { return prob[dir == '+' ? 0 : 1]; }

	double getEstimatedProb(char dir) { return ss[dir == '+' ? 0 : 1] / (ss[0] + ss[1]); }

	double getLogProb(char dir) { return logprob[dir == '+' ? 0 : 1]; }

	//dir must be either + or -
	void update(char dir) { ++ss[dir == '+' ? 0 : 1]; }

	void clear();
	void collect(const Orientation* o);
	void finish(); // Only used if we want to estimate orientation from data

	void read(std::ifstream& fin, int choice); // choice, 0 prob; 1 ss
	void write(std::ofstream& fout, int choice);

	char simulate(Sampler* sampler) { return (sampler->random() < prob[0] ? '+' : '-'); }
	
 private:
 	model_mode_type mode;
 	double *prob, *logprob, *ss;

 	void p2logp();
};

#endif /* ORIENTATION_H_ */
