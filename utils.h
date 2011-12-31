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
#include<cerrno>

const int STRLEN = 10005 ;
const double EPSILON = 1e-300;
const double MINEEL = 1.0;
const double ORIVALVE = 0.1;
const int RANGE = 201;
const int OLEN = 25; // overlap length, number of bases must not be in poly(A) tails
const int NBITS = 32; // use unsigned int, 32 bits per variable

bool verbose = true; // show detail intermediate outputs

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

void printTimeUsed(const time_t& a, const time_t& b, const char* filename = "") {
	int hh = (b - a) / 3600;
	int mm = (b - a) % 3600 / 60;
	int ss = (b - a) % 60;

	printf("Time Used : %d h %02d m %02d s\n", hh, mm, ss);

	if (strcmp(filename, "")) {
		FILE *fo = fopen(filename, "w");
		fprintf(fo, "Time Used : %d h %02d m %02d s\n", hh, mm, ss);
		fclose(fo);
	}
}

void genReadFileNames(const char* readFN, int tagType, int read_type, int& s, char readFs[][STRLEN]){
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

void exitWithError(const char* errmsg) {
	fprintf(stderr, "%s\n", errmsg);
	exit(-1);
}

void pthread_exception(int rc) {
	switch(rc) {
	case EAGAIN:
		fprintf(stderr, "Error code: EAGAIN. Insufficient resources to create another thread, or a system-imposed limit on the number of threads was encountered.\n");
		break;
	case EINVAL:
		fprintf(stderr, "Error code: EINVAL. Invalid settings in attr if pthread_create() is called. Or the implementation has detected that the value specified by thread_id does not refer to a joinable thread if pthread_join() is called.\n");
		break;
	case EPERM:
		fprintf(stderr, "Error code: EPERM. No permission to set the scheduling policy and parameters specified in attr.\n");
		break;
	case EDEADLK:
		fprintf(stderr, "Error code: EDEADLK. A deadlock was detected (e.g., two threads tried to join with each other); or thread_id specifies the calling thread.");
		break;
	case ESRCH:
		fprintf(stderr, "Error code: ESRCH. No thread with thread_id could be found.\n");
		break;
	default: fprintf(stderr, "Unknown error code: %d\n", rc);
	}
	exit(-1);
}

#endif
