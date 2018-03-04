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

#include <cmath>
#include <cassert>
#include <string>
#include <fstream>

#include "MateLenDist.hpp"

MateLenDist::MateLenDist(){
	lb = 1; 
	ub = span = 0;
	ss.clear();
}

void MateLenDist::finish() {
	double sum;

	while (lb <= ub && pmf[lb] == 0.0) ++lb;
	assert(lb <= ub);
	
	sum = 0.0;
	span = ub - lb + 1;
	for (int i = lb; i <= ub; ++i) {
		sum += ss[i];
		ss[i - lb] = ss[i];
	}
	ss.resize(span);
	assert(sum > 0.0);

	pmf = ss;
	for (int i = 0; i < span; ++i) pmf[i] /= sum;
}

double MateLenDist::calcLogP() const {
	double logp = 0.0;
	for (int i = 0; i < span; ++i)
		if (ss[i] > 0.0) logp += ss[i] * log(pmf[i]);
	return logp;
}

void MateLenDist::read(std::ifstream& fin) {
	std::string line;

	assert(fin>> lb>> ub>> span);
	for (int i = 0; i < span; ++i) assert(fin>> pmf[i]);
	getline(fin, line);
}

void MateLenDist::write(std::ofstream& fout, int mate, bool isProb) {
	fout<< "#mld"<< mate<< " 4 MateLenDist, format: lb ub span; [lb, ub], span = ub - lb + 1, probability mass function values"<< std::endl;

	fout<< lb<< '\t'<< ub<< '\t'<< span<< std::endl;  
	for (int i = 1; i < span - 1; ++i) fout<< (isProb ? pmf[i] : ss[i])<< '\t';
	fout<< (isProb ? pmf[span - 1] : ss[span - 1])<< std::endl<< std::endl<< std::endl;
}

void MateLenDist::prepare_for_simulation() {
	cdf = pmf;
	for (int i = 1; i < span; ++i) cdf[i] += cdf[i - 1]; 
}
