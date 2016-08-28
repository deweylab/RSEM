#ifndef UTILS
#define UTILS

#include<cmath>
#include<ctime>
#include<cstdio>
#include<cctype>
#include<cstdlib>
#include<cstring>
#include<cassert>
#include<string>
#include<vector>
#include<stdint.h>

typedef uint64_t HIT_INT_TYPE;
typedef uint64_t READ_INT_TYPE;

const int STRLEN = 10005 ;
const double EPSILON = 1e-300;
const double MINEEL = 1.0;
const double ORIVALVE = 0.1;
const int RANGE = 201;
const int OLEN = 25; // overlap length, number of bases must not be in poly(A) tails
const int NBITS = 32; // use unsigned int, 32 bits per variable

const int MAX_WARNS = 50; // Display at most 50 warnings of the same type

extern bool verbose; // show detail intermediate outputs

inline bool isZero(double a) { return fabs(a) < 1e-8; }
inline bool isLongZero(double a) { return fabs(a) < 1e-30; }

// Assume char's range is -128..127
const int CHAR_RANGE = 128;

static std::vector<int> init_base2id() {
  std::vector<int> vec(CHAR_RANGE, -1);
  vec['a'] = vec['A'] = 0;
  vec['c'] = vec['C'] = 1;
  vec['g'] = vec['G'] = 2;
  vec['t'] = vec['T'] = 3;
  vec['n'] = vec['N'] = 4;

  return vec;
}

static const std::vector<int> base2id = init_base2id();

inline int get_base_id(char c) {
  if (c < 0 || base2id[c] < 0) {
    fprintf(stderr, "Found unknown sequence letter %c at function get_base_id!\n", c);
    exit(-1);
  }
  return base2id[c];
}

static std::vector<int> init_rbase2id() {
  std::vector<int> vec(CHAR_RANGE, -1);
  vec['a'] = vec['A'] = 3;
  vec['c'] = vec['C'] = 2;
  vec['g'] = vec['G'] = 1;
  vec['t'] = vec['T'] = 0;
  vec['n'] = vec['N'] = 4;

  return vec;
}

static const std::vector<int> rbase2id = init_rbase2id();

inline int get_rbase_id(char c) {
  if (c < 0 || rbase2id[c] < 0) {
    fprintf(stderr, "Found unknown sequence letter %c at function get_rbase_id!\n", c);
    exit(-1);
  }
  return rbase2id[c];
}

inline char getOpp(char c) {
  switch(c) {
  case 'a' : return 't';
  case 'c' : return 'g';
  case 'g' : return 'c';
  case 't' : return 'a';
  case 'n' : return 'n';
  case 'A' : return 'T';
  case 'C' : return 'G';
  case 'G' : return 'C';
  case 'T' : return 'A';
  case 'N' : return 'N';
  default :
	fprintf(stderr, "Found unknown sequence letter %c!\n", c);
	exit(-1);
  }
}

inline char getCharacter(int id) {
  switch(id) {
  case 0 : return 'A';
  case 1 : return 'C';
  case 2 : return 'G';
  case 3 : return 'T';
  case 4 : return 'N';
  default :
	  fprintf(stderr, "Found unknown id %d!\n", id);
	  exit(-1);
  }
}

static std::vector<unsigned int> init_mask_code() {
  std::vector<unsigned int> vec(NBITS);
  for (int i = 0; i < NBITS; i++) vec[i] = 1 << i;
  return vec;
}

static std::vector<unsigned int> mask_codes = init_mask_code();

inline std::string cleanStr(const std::string& str) {
  int len = str.length();
  int fr, to;

  fr = 0;
  while (fr < len && isspace(str[fr])) ++fr;
  to = len - 1;
  while (to >= 0 && isspace(str[to])) --to;

  return (fr <= to ? str.substr(fr, to - fr + 1) : "");
}

inline void genReadFileNames(const char* readFN, int tagType, int read_type, int& s, char readFs[][STRLEN]){
	const char tags[3][STRLEN] = {"un", "alignable", "max"};
	char suffix[STRLEN];

	if (read_type == 0 || read_type == 2) {
		strcpy(suffix, "fa");
	}
	else {
		strcpy(suffix, "fq");
	}

	if (read_type == 0 || read_type == 1) {
		s = 1;
		sprintf(readFs[0], "%s_%s.%s", readFN, tags[tagType], suffix);
	}
	else {
		s = 2;
		sprintf(readFs[0], "%s_%s_1.%s", readFN, tags[tagType], suffix);
		sprintf(readFs[1], "%s_%s_2.%s", readFN, tags[tagType], suffix);
	}
}

inline void printTimeUsed(const time_t& a, const time_t& b, const char* program_name) {
	int hh = (b - a) / 3600;
	int mm = (b - a) % 3600 / 60;
	int ss = (b - a) % 60;

	printf("Time Used for %s : %d h %02d m %02d s\n", program_name, hh, mm, ss);
}

inline std::string assemble_command(int argc, char* argv[]) {
  std::string command = argv[0];
  for (int i = 1; i < argc; ++i)
    command += " " + std::string(argv[i]);
  return command;
}

#endif
