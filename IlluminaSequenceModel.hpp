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

#ifndef IlluminaSequenceModel_H_
#define IlluminaSequenceModel_H_

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

#include "MateLenDist.hpp"
#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"
#include "QualDist.hpp"
#include "NoiseProfile.hpp"
#include "NoiseQProfile.hpp"


class IlluminaSequenceModel {
public:
	// mode: 0, master; 1, child; 2 simulation
	// status: 2 bits, bit 1: has qual, bit 2: is first pass
	// min_len: minimum read length; max_len: maximum read length; 
	IlluminaSequenceModel(int mode, int status, int min_len, int max_len); 
	~IlluminaSequenceModel();
	
	/*
		@param   dir             strand information
		@param   pos             position in dir strand, 0-based
		@param   refseq          reference sequence in '+' strand
		@param   cigar           cigar string, as in dir strand 
		@param   seq             read sequence, as in dir strand
		@param   qual            optional quality score sequence, as in dir strand
		@return  log probability of generating the read sequence, no QD
	*/
	double getLogProb(char dir, int pos, const RefSeq& refseq, const CIGARstring& cigar, const SEQstring& seq, const QUALstring& qual) const;

	// for noise reads, no QD
	double getLogProb(const SEQstring& seq, const QUALstring& qual) const {
		return mld->getLobProb(seq.getLen()) + ((status & 1) ? nqpro->getLogProb(qual, seq) : npro->getLogProb(seq));
	}

	/*
		@param   frac            the expected fraction of this alignment
		@param   dir             strand information
		@param   pos             position in dir strand, 0-based
		@param   refseq          reference sequence in '+' strand
		@param   cigar           cigar string
		@param   seq             read sequence
		@param   qual            optional quality score sequence
		@comment   do not update mld qd 
	*/
	void update(double frac, char dir, int pos, const RefSeq& refseq, const CIGARstring& cigar, const SEQstring& seq, const QUALstring& qual);

	// Important, must call for each read, update MLD, QD, NQPRO and NPRO	
	void update(double frac, bool is_aligned, const SEQstring& seq, const QUALstring& qual) {
		if (mode == 0 && (status & 2)) {
			mld->update(seq.getLen());
			if (status & 1) qd->update(qual);
		}
		if (status & 1) nqpro->update(qual, seq, is_aligned, frac);
		else npro->update(seq, is_aligned, frac);
	}

	void clear();
	void collect(const IlluminaSequenceModel* o);
	void finish();
	
	void read(std::ifstream& fin, int choice); // choice: 0 -> p; 1 -> ss; 2 -> c; 
	void write(std::ofstream& fout, int choice);

	// return number of mapped bases, -1 means failed
	int simulate(Sampler* sampler, char dir, int pos, const RefSeq& refseq, std::string& cigar, std::string& readseq, std::string& qual);

private:
	int mode; // 0, master; 1, child; 2, simulation
	int status; // 2 bits, bit 1: has qual, bit 2: is first pass
	int min_len, max_len;

	MateLenDist *mld;
	Markov *markov;
	Profile *pro;
	QProfile *qpro;
	QualDist *qd;
	NoiseProfile *npro;
	NoiseQProfile *nqpro;

	void push_cigar_operation(std::string& cigar, char opchr, int oplen) {
		int s = 0;
		char arr[50];
		
		while (oplen > 0) { arr[s++] = oplen % 10 + '0'; oplen /= 10; }
		while (s > 0) cigar.push_back(arr[--s]);
		cigar.push_back(opchr);
	}
};

// No QualDist prob because it is the same for all alignments of the same read
inline double IlluminaSequenceModel::getLogProb(char dir, int pos, const RefSeq& refseq, const CIGARstring& cigar, const SEQstring& seq, const QUALstring& qual) const {
	double log_prob = mld->getLogProb(seq.getLen());
	int len = cigar.getLen();
	int readpos = 0;

	int prev_state, curr_state; // current Markov state
	char opchr;
	int oplen;

	int ref_base, read_base;

	prev_state = Markov::Begin;
	for (int i = 0; i < len; ++i) {
		opchr = cigar.opchrAt(i);
		oplen = cigar.oplenAt(i);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			// set curr_state
			curr_state = Markov::M;
			// calculate read generating probabilities
			for (int j = 0; j < oplen; ++j) {
				ref_base = refseq.baseCodeAt(dir, pos);
				read_base = seq.baseCodeAt(readpos);
				log_prob += ((status & 1) ? qpro->getLogProb(qual.qualAt(readpos), ref_base, read_base) : pro->getLogProb(readpos, ref_base, read_base));
				++pos; ++readpos;
			}
		}
		else if (opchr == 'I' || opchr == 'S') {
			// set curr_state
			curr_state = (opchr == 'I' ? Markov::I : (prev_state <= Markov::S1 ? Markov::S1 : Markov::S2));
			// calculate read generating probabilities with I/S states
			for (int j = 0; j < oplen; ++j) {
				log_prob += markov->getLogProbBase(curr_state, seq.baseCodeAt(readpos));
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

inline void IlluminaSequenceModel::update(double frac, char dir, int pos, const RefSeq& refseq, const CIGARstring& cigar, const SEQstring& seq, const QUALstring& qual) {
	int len = cigar.getLen();
	int readpos = 0;

	int prev_state, curr_state; // current Markov state
	char opchr;
	int oplen;

	int ref_base, read_base;

	prev_state = Markov::Begin;
	for (int i = 0; i < len; ++i) {
		opchr = cigar.opchrAt(i);
		oplen = cigar.oplenAt(i);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			// set curr_state
			curr_state = Markov::M;
			// update read generating probabilities
			for (int j = 0; j < oplen; ++j) {
				ref_base = refseq.baseCodeAt(dir, pos);
				read_base = seq.baseCodeAt(readpos);
				if (status & 1) qpro->update(qual.qualAt(readpos), ref_base, read_base, frac);
				else pro->update(readpos, ref_base, read_base, frac);
				++pos; ++readpos;
			}
		}
		else if (opchr == 'I' || opchr == 'S') {
			// set curr_state
			curr_state = (opchr == 'I' ? Markov::I : (prev_state <= Markov::S1 ? Markov::S1 : Markov::S2));
			// update read generating probabilities with I/S states
			for (int j = 0; j < oplen; ++j) {
				markov->updateBase(curr_state, seq.baseCodeAt(readpos), frac);
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

// assume pos is valid
inline int IlluminaSequenceModel::simulate(Sampler* sampler, char dir, int pos, const RefSeq& refseq, std::string& cigar, std::string& readseq, std::string& qual) {
	int curr_state; // current Markov states
	int qval, refpos, readpos, reflen, readlen;
	
	char opchr, curr_opchr, non_s_opchr;
	int oplen;
	
	qval = -1; 
	refpos = pos; readpos = 0;
	reflen = refseq.getLen();
	readlen = mld->simulate(sampler); // simulate read length
	curr_opchr = 0; oplen = 0; non_s_opchr = 0;
	cigar = readseq = qual = "";

	curr_state = Markov::Begin;
	do {
		curr_state = markov->simulate(sampler, curr_state);
		opchr = Markov::getOpChr(curr_state);
		if (opchr != 'S') non_s_opchr = opchr;

		if (curr_opchr != opchr) {
			if (oplen > 0) push_cigar_operation(cigar, curr_opchr, oplen);
			curr_opchr = opchr;
			oplen = 0;
		}
		++oplen;
		
		if (opchr == 'M') {
			if (refpos >= reflen) return -1; // failed

			if (stats & 1) {
				qval = (qval < 0 ? qd->simulate(sampler) : qd->simulate(sampler, qval));
				qual.push_back(qval2char(qval));
				readseq.push_back(qpro->simulate(sampler, qval, refseq.baseCodeAt(dir, refpos)));
			}
			else readseq.push_back(pro->simulate(sampler, readpos, refseq.baseCodeAt(dir, refpos)));

			++refpos; ++readpos;
		}
		else if (opchr == 'S' || opchr == 'I') {
			if (status & 1) {
				qval = (qval < 0 ? qd->simulate(sampler) : qd->simulate(sampler, qval));
				qual.push_back(qval2char(qval));
			}
			readseq.push_back(markov->simulateBase(sampler, curr_state));
			++readpos;
		}
		else {
			assert(opchr == 'D');
			if (refpos >= reflen) return -1;
			++refpos;
		}
	} while (readpos < readlen); 

	if (non_s_opchr != 'M') return -1;

	if (oplen > 0) push_cigar_operation(cigar, curr_opchr, oplen);

	return refpos - pos;
}

#endif /*IlluminaSequenceModel*/
