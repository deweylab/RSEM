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

#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include <cctype>
#include <cstdint>
#include <string>
#include <vector>

typedef uint64_t HIT_INT_TYPE;
typedef uint64_t READ_INT_TYPE;

extern bool verbose; // show detail intermediate outputs

const std::string VERSION = "2.0";

const int STRLEN = 10005 ;
const double EPSILON = 1e-300;
const int MAX_WARNS = 50; // Display at most 50 warnings of the same type

const double MINEEL = 1.0;
const double ORIVALVE = 0.1;
const int RANGE = 201;
const int OLEN = 25; // overlap length, number of bases must not be in poly(A) tails
const int NBITS = 32; // use unsigned int, 32 bits per variable
const int MASK_LEN = 24; // the last MASK_LEN bp of a sequence cannot be aligned if poly(A) tail is added

const int NCODES = 5; // A, C, G, T, N
const int QSIZE = 100; // quality score range, from 0 to QSIZE - 1

inline char qval2char(int qval) { return char(qval + 33); }
inline int char2qval(int c) { return c - 33; }  

/*
	In our context, isXXZero is called only for non-negative values. Thus we omit the fabs function
 */
//inline bool isZero(double a) { return fabs(a) < 1e-8; }
//inline bool isLongZero(double a) { return fabs(a) < 1e-30; }
inline bool isZero(double a) { return a < 1e-8; }
inline bool isLongZero(double a) { return a < 1e-30; }

inline std::string cleanStr(const std::string& str) {
	int len = str.length();
	int fr, to;

	fr = 0;
	while (fr < len && isspace(str[fr])) ++fr;
	to = len - 1;
	while (to >= 0 && isspace(str[to])) --to;

	return (fr <= to ? str.substr(fr, to - fr + 1) : "");
}

inline std::string generateCommand(int argc, char* argv[]) {
	std::string command = "\"" + std::string(argv[0]);
	for (int i = 1; i < argc; ++i) command += " " + std::string(argv[i]);
	command += "\"";

	return command;
}



static const double pseudoC[NCODES] = {1.0, 1.0, 1.0, 1.0, 0.1}; // pseudo count for A/C/G/T/N
static const double initP[4] = {0.9, 0.03, 0.01, 0.2}; // 0, ref_base == read_base; 1, ref_base != read_base & read_base != 'N'; 2, read_base == 'N'; 3, ref_base == 'N', uniform

static const char code2base[] = "ACGTN";

static std::vector<char> init_base2rbase() {
	std::vector<char> vec(128, -1);
	vec['a'] = 't'; vec['A'] = 'T';
	vec['c'] = 'g'; vec['C'] = 'G';
	vec['g'] = 'c'; vec['G'] = 'C';
	vec['t'] = 'a'; vec['T'] = 'A';
	vec['n'] = 'n'; vec['N'] = 'N';

	return vec;
}

static const std::vector<char> base2rbase = init_base2rbase();

static std::vector<int> init_base2code() {
	std::vector<int> vec(128, -1);
	vec['a'] = vec['A'] = 0;
	vec['c'] = vec['C'] = 1;
	vec['g'] = vec['G'] = 2;
	vec['t'] = vec['T'] = 3;
	vec['n'] = vec['N'] = 4;
	
	return vec;
}

static const std::vector<int> base2code = init_base2code();

static std::vector<int> init_rbase2code() {
	std::vector<int> vec(128, -1);
	vec['a'] = vec['A'] = 3;
	vec['c'] = vec['C'] = 2;
	vec['g'] = vec['G'] = 1;
	vec['t'] = vec['T'] = 0;
	vec['n'] = vec['N'] = 4;
	
	return vec;
}

static const std::vector<int> rbase2code = init_rbase2code();

#endif
