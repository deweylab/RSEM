#ifndef __GTFITEM__
#define __GTFITEM__

#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<string>
#include<sstream>

#include "utils.h"

class GTFItem {
public:

	GTFItem() {
		seqname = source = feature = "";
		score = "";
		start = end = 0;
		strand = 0; //strand is a char variable
		frame = "";
		gene_id = transcript_id = "";
		left = "";
	}

	bool operator<(const GTFItem& o) const {
		if (gene_id != o.gene_id) return gene_id < o.gene_id;
		if (transcript_id != o.transcript_id) return transcript_id < o.transcript_id;
		return start < o.start;
	}

	void my_assert(char value, std::string& line, const std::string& msg) {
		if (!value) {
			fprintf(stderr, ".gtf file might be corrupted!\n");
			fprintf(stderr, "Stop at line : %s\n", line.c_str());
			fprintf(stderr, "Error Message: %s\n", msg.c_str());
			exit(-1);
		}
	}

	void parse(std::string line) {
		std::istringstream strin(line);
		std::string tmp;

		getline(strin, seqname, '\t');
		getline(strin, source, '\t');
		getline(strin, feature, '\t');
		getline(strin, tmp, '\t');
		start = atoi(tmp.c_str());
		getline(strin, tmp, '\t');
		end = atoi(tmp.c_str());
		getline(strin, score, '\t');
		getline(strin, tmp, '\t');
		my_assert((tmp.length() == 1 && (tmp[0] == '+' || tmp[0] == '-')), line, "Strand is neither '+' nor '-'!");
		strand = tmp[0];
		getline(strin, frame, '\t');

		getline(strin, left); // assign attributes and possible comments into "left"

		strin.clear(); strin.str(left);
		bool find_gene_id = false, find_transcript_id = false;

		while (getline(strin, tmp, ';') && (!find_gene_id || !find_transcript_id)) {
			tmp = cleanStr(tmp);
			size_t pos = tmp.find(' ');
			my_assert((pos != std::string::npos), line, "Cannot separate the identifier from the value for attribute " + tmp + "!");
			std::string identifier = tmp.substr(0, pos);

			if (identifier == "gene_id") {
				my_assert(!find_gene_id, line, "gene_id appear more than once!");
				tmp = cleanStr(tmp.substr(pos));
				my_assert((tmp[0] == '"' && tmp[tmp.length() - 1] == '"'), line, "Textual attributes should be surrounded by doublequotes!");
				gene_id = tmp.substr(1, tmp.length() - 2);
				find_gene_id = true;
			} else if (identifier == "transcript_id") {
				my_assert(!find_transcript_id, line, "transcript_id appear more than once!");
				tmp = cleanStr(tmp.substr(pos));
				my_assert((tmp[0] == '"' && tmp[tmp.length() - 1] == '"'), line, "Textual attributes should be surrounded by doublequotes!");
				transcript_id = tmp.substr(1, tmp.length() - 2);
				find_transcript_id = true;
			}
		}

		my_assert(feature != "exon" || find_gene_id, line, "Cannot find gene_id!");
		my_assert(feature != "exon" || find_transcript_id, line, "Cannot find transcript_id!");
		if (!find_gene_id && feature != "exon") { printf("Warning: line \" %s \" does not contain a gene_id attribute! Since this line will not be used for reference construction, it is skipped. But if you think this GTF file is corrupted, you should find a complelete GTF file instead and rebuild the reference.\n", line.c_str()); }
		if (!find_transcript_id && feature != "exon") { printf("Warning: line \" %s \" does not contain a transcript_id attribute! Since this line will not be used for reference construction, it is skipped. But if you think this GTF file is corrupted, you should find a complelete GTF file instead and rebuild the reference.\n", line.c_str()); }
	}

	std::string getSeqName() { return seqname; }
	std::string getSource() { return source; }
	std::string getFeature() { return feature; }
	int getStart() { return start; }
	int getEnd() { return end; }
	char getStrand() { return strand; }
	std::string getScore() { return score; }  // float, integer or "." ; let downstream programs parse it
	std::string getFrame() { return frame; }  // 0, 1, 2, or "."; let downstream programs parse it
	std::string getGeneID() { return gene_id; }
	std::string getTranscriptID() { return transcript_id; }
	std::string getLeft() { return left; }

	void setGeneID(const std::string& gene_id) {
		this->gene_id = gene_id;
	}

	std::string toString() {
		std::string val;
		std::ostringstream strout;
		strout<<seqname<<'\t'<<source<<'\t'<<feature<<'\t'<<start<<'\t'<<end<<'\t'<<score<<'\t'<<strand<<'\t'<<frame<<'\t';
		strout<<"gene_id \""<<gene_id<<"\"; transcript_id \""<<transcript_id<<"\";"<<left;
		val = strout.str();

		return val;
	}

private:
	std::string seqname, source, feature;
	std::string score;
	int start, end;
	char strand;
	std::string frame;
	std::string gene_id, transcript_id;
	std::string left;
};

#endif
