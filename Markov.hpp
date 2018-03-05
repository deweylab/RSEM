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

#ifndef MARKOV_H_
#define MARKOV_H_

#include <vector>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Markov {
public:
	static const int NSTATES = 6;
	
	static const int Begin = 0; // Begin -> S1, M
	static const int S1 = 1; // S1 -> S1, M
	static const int M = 2; // M -> M, I, D, S2, & End only if M is at the end of the fragment and no S2 observed
	static const int I = 3; // I -> I, M, D
	static const int D = 4;	// D -> M, I, D, & End only if D is at the end of the fragment; for simulation, allow D -> S2 if S2 observed
	static const int S2 = 5;
	static const int End = 6; // This state is defined but not in P
	
	static char getOpChr(int state) { return state2char[state]; }
	
	Markov(bool to_log_space = false);
	~Markov();

	// probabilities are in log space
	double getLogProb(int curr_state, int prev_state = Begin) const { return P[prev_state][curr_state]; }

	// probabilities in log space
	double getLogProbBase(int state, int code) const {
		assert(s2pos[state] >= 0 && code >= 0);
		return probB[s2pos[state]][code];
	}

	void update(int curr_state, double frac, int prev_state = Begin) {
		P[prev_state][curr_state] += frac;
	}

	void updateBase(int state, int code, double frac) {
		assert(s2pos[state] >= 0 && code >= 0);
		probB[s2pos[state]][code] += frac;
	}
	
	void collect(const Markov* o);
	void finish();
	void clear();
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, bool isProb = true);
	
	void prepare_for_simulation();
	
	// only call after prepare_for_simulation is called
	bool hasS2() const { return P[S2][0] > 0.0; }

	int simulate(Sampler *sampler, int state = Begin) {
		return offset[state] + sampler->sample(P[state], lens[state]);
	}
	
	char simulateBase(Sampler *sampler, int state) {
		return code2base[sampler->sample(probB[s2pos[state]], NCODES)];
	}
	
private:
	static const char state2char[NSTATES]; // state to char
	static const int lens[NSTATES]; // number of transition states for each state
	static const int offset[NSTATES]; // offsets	
	static const int s2pos[NSTATES]; // state to row position in probB
	
	double P[NSTATES][NSTATES], (*ss_P)[NSTATES];
	double probB[3][NCODES], (*ss_probB)[NCODES]; // probability of generating bases given state S1, I/I2, S2

	void ss2p();
};

#endif
