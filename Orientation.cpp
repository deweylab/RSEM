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

#include <new>
#include <cmath>
#include <limits>
#include <cassert>
#include <string>
#include <fstream>

#include "Orientation.hpp"

Orientation::Orientation(int mode) : mode(mode), prob(NULL), logprob(NULL), ss(NULL) {
	if (mode != 1) prob = new double[2];
	if (mode == 0) logprob = new double[2];
	if (mode != 2) ss = new double[2];
}

Orientation::~Orientation() {
	if (prob != NULL) delete[] prob;
	if (logprob != NULL) delete[] logprob;
	if (ss != NULL) delete[] ss;
}

void Orientation::clear() { 
	ss[0] = ss[1] = 0.0;
}

void Orientation::collect(const Orientation* o) {
	ss[0] += o->ss[0];
	ss[1] += o->ss[1];
}

void Orientation::finish() {
	double sum = ss[0] + ss[1];
	prob[0] = ss[0] / sum;
	prob[1] = ss[1] / sum;
	p2logp();
}

void read(std::ifstream& fin, int choice) {
	std::string line;
	double *in = NULL;

	switch(choice) {
		case 0: in = prob; break;
		case 1: in = ss;
	}

	assert((fin>> line) && (line == "#ori"));
	getline(fin, line);
	assert(fin>> in[0]>> in[1]);

	if (mode == 0 && choice == 0) p2logp();
}

void write(std::ofstream& fout, int choice) {
	double *out = NULL;

	switch(choice) {
		case 0: out = prob; break;
		case 1: out = ss;
	}

	fout<< "#ori\tOrientation\t3\tformat: forward_probability reverse_probability"<< std::endl;
	fout<< out[0]<< ' '<< out[1]<< std::endl<< std::endl<< std::endl;
}

void Orientation::p2logp() {
	logprob[0] = (prob[0] > 0.0 ? log(prob[0]) : -std::numeric_limits<double>::infinity());
	logprob[1] = (prob[1] > 0.0 ? log(prob[1]) : -std::numeric_limits<double>::infinity());
}
