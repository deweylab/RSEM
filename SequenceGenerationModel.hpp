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

#ifndef SEQUENCINGERRORMODEL_H_
#define SEQUENCINGERRORMODEL_H_

#include <cassert>
#include <string>
#include <ostream>
#include <fstream>
#include <sstream>

#include "utils.h"
#include "sampling.hpp"

#include "RefSeq.hpp"
#include "CIGARstring.hpp"
#include "QUALstring.hpp"
#include "SEQstring.hpp"

#include "Markov.hpp"
#include "Profile.hpp"
#include "QualDist.hpp"
#include "QProfile.hpp"

// This class is an abstraction of the NGS sequencing error model, it is applicable to different types of sequencing, e.g. RNA-Seq, ChIP-Seq, Structure-Seq  
class SequencingErrorModel {
public:
	/*
		@param   hasQual      if we have quality scores
		@param   to_log_space if we should convert probability to log space
		@param   maxL         maximum read/mate length 
	*/
	SequencingErrorModel(bool hasQual, bool to_log_space = false, int maxL = 0);  
	~SequencingErrorModel();
	
	/*
		@param   dir             strand information
		@param   pos             position in dir strand, 0-based
		@param   refseq          reference sequence in '+' strand
		@param   cigar           cigar string
		@param   seq             read sequence
		@param   qual            optional quality score sequence
		@return  log probability of generating the read sequence
	*/
	double getLogProb(char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual = NULL) const;

	/*
		@param   frac            the expected fraction of this alignment
		@param   dir             strand information
		@param   pos             position in dir strand, 0-based
		@param   refseq          reference sequence in '+' strand
		@param   cigar           cigar string
		@param   seq             read sequence
		@param   qual            optional quality score sequence
	*/
	void update(double frac, char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual = NULL);
	
	void collect(const SequencingErrorModel* o);
	void finish(int length = -1); // for parseAlignment.cpp, length set the observed maximum read/mate length
	void clear();
	
	void read(std::ifstream& fin);
	void write(std::ofstream& fout, bool isProb = true);

	void prepare_for_simulation(QualDist* qd);

	void simulate(Sampler* sampler, int frag_len, int len, char dir, int pos, const RefSeq* refseq, std::string& cigar, std::string& readseq, std::string& qual);

private:
	bool hasQual; // if we have quality scores

	Markov *markov;
	Profile *profile;
	QProfile *qprofile;
	
	QualDist *qd; // qd is only used for simulation

	void push_cigar_operation(std::string& cigar, char opchr, int oplen) {
		int s = 0;
		char arr[50];
		
		while (oplen > 0) { arr[s++] = oplen % 10 + '0'; oplen /= 10; }
		while (s > 0) cigar.push_back(arr[--s]);
		cigar.push_back(opchr);
	}
};

// No QualDist prob because it is the same for all alignments of the same read
inline double SequencingErrorModel::getLogProb(char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) const {
	double log_prob = 0.0;
	int len = cigar->getLen();
	int readpos = 0;

	int prev_state, curr_state; // current Markov state
	char opchr;
	int oplen;

	int ref_base, read_base;

	prev_state = Markov::Begin;
	for (int i = 0; i < len; ++i) {
		opchr = cigar->opchrAt(i);
		oplen = cigar->oplenAt(i);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			// set curr_state
			curr_state = Markov::M;
			// calculate read generating probabilities
			for (int j = 0; j < oplen; ++j) {
				ref_base = refseq->baseCodeAt(dir, pos);
				read_base = seq->baseCodeAt(readpos);
				log_prob += (hasQual ? qprofile->getLogProb(qual->qualAt(readpos), ref_base, read_base) : profile->getLogProb(readpos, ref_base, read_base));
				++pos; ++readpos;
			}
		}
		else if (opchr == 'I' || opchr == 'S') {
			// set curr_state
			curr_state = (opchr == 'I' ? Markov::I : (prev_state <= Markov::S1 ? Markov::S1 : Markov::S2));
			// calculate read generating probabilities with I/S states
			for (int j = 0; j < oplen; ++j) {
				log_prob += markov->getLogProbBase(curr_state, seq->baseCodeAt(readpos));
				++readpos;
			}
		}
		else {
			assert(opchr == 'D');
			// set curr_state
			curr_state = Markov::D;
			pos += oplen;
		}

		// adding up the log probabilities from Markov chain itself
		log_prob += markov->getLogProb(curr_state, prev_state);
		if (oplen > 1) log_prob += markov->getLogProb(curr_state, curr_state) * (oplen - 1);

		prev_state = curr_state;
	}
	
	return log_prob;
}

inline void SequencingErrorModel::update(double frac, char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) {
	int len = cigar->getLen();
	int readpos = 0;

	int prev_state, curr_state; // current Markov state
	char opchr;
	int oplen;

	int ref_base, read_base;

	prev_state = Markov::Begin;
	for (int i = 0; i < len; ++i) {
		opchr = cigar->opchrAt(i);
		oplen = cigar->oplenAt(i);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			// set curr_state
			curr_state = Markov::M;
			// update read generating probabilities
			for (int j = 0; j < oplen; ++j) {
				ref_base = refseq->baseCodeAt(dir, pos);
				read_base = seq->baseCodeAt(readpos);
				if (hasQual) qprofile->update(qual->qualAt(readpos), ref_base, read_base, frac);
				else profile->update(readpos, ref_base, read_base, frac);
				++pos; ++readpos;
			}
		}
		else if (opchr == 'I' || opchr == 'S') {
			// set curr_state
			curr_state = (opchr == 'I' ? Markov::I : (prev_state <= Markov::S1 ? Markov::S1 : Markov::S2));
			// update read generating probabilities with I/S states
			for (int j = 0; j < oplen; ++j) {
				markov->update(curr_state, seq->baseCodeAt(readpos), frac);
				++readpos;
			}
		}
		else {
			assert(opchr == 'D');
			// set curr_state
			curr_state = Markov::D;
			pos += oplen;
		}

		// update Markov chain parameters
		markov->update(curr_state, frac, prev_state);
		if (oplen > 1) markov->update(curr_state, frac * (oplen - 1), curr_state);

		prev_state = curr_state;
	}
}

inline void SequencingErrorModel::simulate(Sampler* sampler, int fraglen, int readlen, char dir, int pos, const RefSeq* refseq, std::string& cigar, std::string& readseq, std::string& qual) {
	int curr_state; // current Markov states
	int qval = -1;
	
	char opchr, curr_opchr = 0;
	int oplen = 0;
	
	cigar = readseq = qual = "";
	curr_state = Markov::Begin;
	do {
		curr_state = markov->simulate(sampler, curr_state);
		opchr = Markov::getOpChr(curr_state);
		
		if (curr_opchr != opchr) {
			if (oplen > 0) push_cigar_operation(cigar, curr_opchr, oplen);
			curr_opchr = opchr;
			oplen = 0;
		}
		++oplen;
		
		if (opchr == 'M') {
			if (hasQual) {
				qval = (qval < 0 ? qd->simulate(sampler) : qd->simulate(sampler, qval));
				qual.push_back(qval2char(qval));
				readseq.push_back(qprofile->simulate(sampler, qval, refseq->baseCodeAt(dir, pos)));
			}
			else readseq.push_back(profile->simulate(sampler, readpos, refseq->baseCodeAt(dir, pos)));
			--readlen; --fraglen;
		}
		else if (opchr == 'S' || opchr == 'I') {
			if (hasQual) {
				qval = (qval < 0 ? qd->simulate(sampler) : qd->simulate(sampler, qval));
				qual.push_back(qval2char(qval));
			}
			readseq.push_back(markov->simulateBase(sampler, curr_state));
			--readlen;
		}
		else {
			assert(opchr == 'D');
			--fraglen;
		}
	} while (readlen > 0 && fraglen > 0); 

	if (oplen > 0) push_cigar_operation(cigar, curr_opchr, oplen);

	// If has S2, the rest sequence are generated as S2 state
	if (readlen > 0 && markov->hasS2()) {
		push_cigar_operation(cigar, 'S', readlen);
		for (int i = 0; i < readlen; ++i) {
			if (hasQual) {
				qval = (qval < 0 ? qd->simulate(sampler) : qd->simulate(sampler, qval));
				qual.push_back(qval2char(qval));
			}
			readseq.push_back(markov->simulateBase(sampler, Markov::S2));
		}
	}	
}

#endif
