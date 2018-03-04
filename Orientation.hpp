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

#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "sampling.hpp"

class Orientation {
 public:
	Orientation(double probF = 0.5) {
		prob[0] = probF;
		prob[1] = 1.0 - probF;

		ss[0] = ss[1] = 0.0;
	}
	
	//dir : +/-
	double getProb(char dir) { return prob[dir == '+' ? 0 : 1]; }
	
	//dir must be either + or -
	void update(char dir) { ++ss[dir == '+' ? 0 : 1]; }

	void collect(const Orientation* o) {
		ss[0] += o->ss[0];
		ss[1] += o->ss[1];
	}

	// Only used if we want to estimate it from data
	void finish() {
		double sum = ss[0] + ss[1];
		prob[0] = ss[0] / sum;
		prob[1] = ss[1] / sum;
	}

	void clear() { ss[0] = ss[1] = 0.0; }

	void read(std::ifstream& fin) {
		assert(fin>> prob[0]); 
		getline(fin, line);
		prob[1] = 1.0 - prob[0];
	}
	
	void write(std::ofstream& fout, bool isProb) {
		if (isProb) {
			fout<< "#ori 3 Orientation, format: forward probability"<< std::endl;
			fout<< prob[0]<< std::endl;
		}
		else {
			fout<< "#ori 3 Orientation, format: ss[0], ss[1], estimated_forward_probability"<< std::endl;
			fout<< ss[0]<< ss[1]<< ss[0] / (ss[0] + ss[1])<< std::endl;
		}

		fout<< std::endl<< std::endl;
	}
	
	char simulate(Sampler* sampler) { return (sampler->random() < prob[0] ? '+' : '-'); }
	
 private:  
	double prob[2]; // 0 + 1 -
	double ss[2]; // sufficient statistics
};

#endif /* ORIENTATION_H_ */
