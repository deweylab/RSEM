#ifndef REFSEQ
#define REFSEQ

#include<cassert>
#include<fstream>
#include<string>
#include<vector>

#include "utils.h"

//Each Object can only be used once
class RefSeq {
public:
	RefSeq() {
		fullLen = totLen = 0;
		name = ""; seq = "";
		fmasks.clear();
	}

	//Constructor , seq : the forward strand of the reference
	//tag does not contain ">"
	//polyALen : length of polyA tail we add
	RefSeq(const std::string& name, const std::string& seq, int polyALen) {
		fullLen = seq.length();
		totLen = fullLen + polyALen;
		this->name = name;
		this->seq = seq;
		this->seq.append(polyALen, 'A');

		assert(fullLen > 0 && totLen >= fullLen);

		int len = (fullLen - 1) / NBITS + 1;
		fmasks.assign(len, 0);
		// set mask if poly(A) tail is added
		if (polyALen > 0) {
			for (int i = std::max(fullLen - OLEN + 1, 0); i < fullLen; i++) setMask(i);
		}
  }

	RefSeq(const RefSeq& o) {
		fullLen = o.fullLen;
		totLen = o.totLen;
		name = o.name;
		seq = o.seq;
		fmasks = o.fmasks;
	}

	RefSeq& operator= (const RefSeq &rhs) {
		if (this != &rhs) {
			fullLen = rhs.fullLen;
			totLen = rhs.totLen;
			name = rhs.name;
			seq = rhs.seq;
			fmasks = rhs.fmasks;
		}

		return *this;
	}

	~RefSeq() {}

	bool read(std::ifstream&, int  = 0);
	void write(std::ofstream&);

	int getFullLen() const { return fullLen; }

	int getTotLen() const { return totLen; }

	const std::string& getName() const { return name; }

	std::string getSeq() const { return seq; }

	std::string getRSeq() const {
		std::string rseq = "";
		for (int i = totLen - 1; i >= 0; i--) rseq.push_back(getCharacter(get_rbase_id(seq[i])));
		return rseq;
	}

	//get the sequence  dir 0 : + 1 : -
	std::string getSeq(int dir) const {
		return (dir == 0 ? getSeq() : getRSeq());
	}
  
	int get_id(int pos, int dir) const {
		assert(pos >= 0 && pos < totLen);
		return (dir == 0 ? get_base_id(seq[pos]) : get_rbase_id(seq[totLen - pos - 1]));
	}

	bool getMask(int seedPos) const {
		assert(seedPos >= 0 && seedPos < totLen);
		return fmasks[seedPos / NBITS] & mask_codes[seedPos % NBITS];
	}

	void setMask(int seedPos) {
		assert(seedPos >= 0 && seedPos < totLen);
		fmasks[seedPos / NBITS] |= mask_codes[seedPos % NBITS];
	}
  
private:
	int fullLen; // fullLen : the original length of an isoform
	int totLen; // totLen : the total length, included polyA tails, if any
	std::string name; // the tag
	std::string seq; // the raw sequence, in forward strand
	std::vector<unsigned int> fmasks; // record masks for forward strand, each position occupies 1 bit
};

//internal read; option 0 : read all 1 : do not read seqences
bool RefSeq::read(std::ifstream& fin, int option) {
	std::string line;

	if (!(fin>>fullLen>>totLen)) return false;
	assert(fullLen > 0 && totLen >= fullLen);
	getline(fin, line);
	if (!getline(fin, name)) return false;
	if (!getline(fin, seq)) return false;

	int len = (fullLen - 1) / NBITS + 1; // assume each cell contains NBITS bits
	fmasks.assign(len, 0);
	for (int i = 0; i < len; i++)
	    if (!(fin>>fmasks[i])) return false;
	getline(fin, line);

	assert(option == 0 || option == 1);
	if (option == 1) { seq = ""; }

	return true;
}

//write to file in "internal" format
void RefSeq::write(std::ofstream& fout) {
	fout<<fullLen<<" "<<totLen<<std::endl;
	fout<<name<<std::endl;
	fout<<seq<<std::endl;

	int len = fmasks.size();
	for (int i = 0; i < len - 1; i++) fout<<fmasks[i]<<" ";
	fout<<fmasks[len - 1]<<std::endl;
}

#endif
