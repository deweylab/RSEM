#ifndef POLYARULES
#define POLYARULES

#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<set>
#include<cstring>
#include<fstream>

/**
Isoform id starts from 1 !
*/

class PolyARules {
public:
	PolyARules() {
		polyAChoice = 0;
		polyALen = 0;
		exceptionList.clear();
	}

	//Assume parameters are valid here
	PolyARules(int polyAChoice, int polyALen, char* exceptionF) {
		this->polyAChoice = polyAChoice;
		this->polyALen =  polyALen;

		if (polyAChoice == 2) {
			exceptionList.clear();

			std::string transcript_id;
			std::ifstream fin(exceptionF);
			if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", exceptionF); exit(-1); }

			while (fin>> transcript_id) {
				exceptionList.insert(transcript_id);
			}

			fin.close();
		}
	}

	//get the length of padding poly As
	int getLenAt(const std::string& transcript_id) {
		switch(polyAChoice) {
		case 0 : return polyALen;
		case 1 : return 0;
		case 2 : iter = exceptionList.find(transcript_id);
				 return (iter == exceptionList.end() ? polyALen : 0);
		default : assert(false);
		}
	}

private:
	int polyAChoice; // 0, pad; 1, do not pad; 2 pad all but those in exceptionList
	int polyALen;
	std::set<std::string> exceptionList; // exception list of transcript_ids
	std::set<std::string>::iterator iter;
};

#endif
