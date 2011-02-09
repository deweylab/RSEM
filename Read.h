#ifndef READ
#define READ

/**
father class of SingleRead, SingleReadQ, PairedEndRead, PairedEndReadQ
 */

#include<iostream>
#include<string>

class Read {
	public:
		Read() { name = ""; low_quality = false; }
		bool read(int argc, std::istream* argv[], int flags = 7) { return false; }  //read from file, flags, which entries loaded 1 : readseq, 2 : quality score 4 : name
		void write(int argc, std::ostream* argv[]) {}; //write to files // do not write if does not read fully
		const std::string& getName() const { return name; }
		bool isLowQuality() const { return low_quality; } // if this read is low quality and should not be used
	protected:
		std::string name; //name of the read
		bool low_quality;
};

#endif
