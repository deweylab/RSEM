/* Copyright (c) 2015
	 Bo Li (University of California, Berkeley)
	 bli25@berkeley.edu

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

#ifndef GTFITEM_H_
#define GTFITEM_H_

#include<cstdio>
#include<cctype>
#include<cstdlib>
#include<cassert>
#include<string>
#include<sstream>

class GTFItem {
 public:
	
	GTFItem() {
		seqname = source = feature = "";
		score = "";
		start = end = 0;
		strand = 0; //strand is a char variable
		frame = "";
		gene_id = transcript_id = "";
		gene_name = transcript_name = "";
		left = "";
	}
	
	bool operator<(const GTFItem& o) const {
		if (gene_id != o.gene_id) return gene_id < o.gene_id;
		if (transcript_id != o.transcript_id) return transcript_id < o.transcript_id;
		return start < o.start;
	}
	
	void parse(const std::string& line) {
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
		gtf_assert((tmp.length() == 1 && (tmp[0] == '+' || tmp[0] == '-')), line, "Strand is neither '+' nor '-'!");
		strand = tmp[0];
		getline(strin, frame, '\t');
		
		getline(strin, left); // assign attributes and possible comments into "left"
	}

	void parseAttributes(const std::string& line) {
		assert(feature == "exon");
		gene_id = transcript_id = "";
		gene_name = transcript_name = "";

		int nleft = 4;
		int pos, lpos = 0, rpos, left_len = left.length();    
		std::string identifier;

		while (nleft > 0 && get_an_attribute(lpos, rpos, left_len)) {
			// locate identifier
			pos = lpos;
			while (pos < rpos && !isspace(left[pos])) ++pos;
			gtf_assert(isspace(left[pos]), line, "Cannot locate the identifier from attribute " + left.substr(lpos, rpos + 1 - lpos) + "!");
			identifier = left.substr(lpos, pos - lpos);

			// prepare for the next attribute
			lpos = rpos + 1;

			// locate value
			++pos;
			while (pos < rpos && isspace(left[pos])) ++pos;
			if (left[pos] != '"') pos = rpos;
			--rpos;
			while (rpos > pos && isspace(left[rpos])) --rpos;
			if (rpos > pos && left[rpos] != '"') rpos = pos;

			// test if the identifier is interested
			if (identifier == "gene_id") {
				gtf_assert(gene_id == "", line, "gene_id appear more than once!");
				gtf_assert(rpos - pos > 1 , line, "Attribute " + identifier + "'s value should be surrounded by double quotes and cannot be empty!");
				gene_id = left.substr(pos + 1, rpos - pos - 1);
				--nleft;
			}
			else if (identifier == "transcript_id") {
				gtf_assert(transcript_id == "", line, "transcript_id appear more than once!");
				gtf_assert(rpos - pos > 1 , line, "Attribute " + identifier + "'s value should be surrounded by double quotes and cannot be empty!");
				transcript_id = left.substr(pos + 1, rpos - pos - 1);
				--nleft;
			}
			else if (identifier == "gene_name" && gene_name == "" && rpos - pos > 1) {        
				gene_name = left.substr(pos + 1, rpos - pos - 1);
				--nleft;
			}
			else if (identifier == "transcript_name" && transcript_name == "" && rpos - pos > 1) {
				transcript_name = left.substr(pos + 1, rpos - pos - 1);
				--nleft;
			}
		}
		
		gtf_assert(gene_id != "", line, "Cannot find gene_id!");
		gtf_assert(transcript_id != "", line, "Cannot find transcript_id!");
	}
	
	const std::string& getSeqName() const { return seqname; }
	const std::string& getSource() const { return source; }
	const std::string getFeature() const { return feature; }
	int getStart() const { return start; }
	int getEnd() const { return end; }
	char getStrand() const { return strand; }
	const std::string& getScore() const { return score; }  // float, integer or "." ; let downstream programs parse it
	const std::string& getFrame() const { return frame; }  // 0, 1, 2, or "."; let downstream programs parse it
	const std::string& getGeneID() const { return gene_id; }
	const std::string& getTranscriptID() const { return transcript_id; }
	const std::string& getGeneName() const { return gene_name; }
	const std::string& getTranscriptName() const { return transcript_name; }
	const std::string getLeft() { return left; }
	
	void setGeneID(const std::string& gene_id) {
		this->gene_id = gene_id;
	}
	
	std::string toString() {
		std::ostringstream strout("");
		strout<< seqname<< '\t'<< source<< '\t'<< feature<< '\t'<< start<< '\t'<< end<< '\t'<< score<< '\t'<< strand<< '\t'<< frame<< '\t'<< left;
		return strout.str();
	}
	
 private:
	std::string seqname, source, feature;
	std::string score;
	int start, end;
	char strand;
	std::string frame;
	std::string gene_id, transcript_id;
	std::string gene_name, transcript_name;
	std::string left;


	bool get_an_attribute(int& lpos, int& rpos, int left_len) {
		bool in_quote;

		while (lpos < left_len && isspace(left[lpos])) ++lpos; // remove leading white spaces

		rpos = lpos; in_quote = false;
		while (rpos < left_len && (left[rpos] != ';' || in_quote)) {
			if (left[rpos] == '"') in_quote ^= true;
			++rpos;
		}

		return rpos < left_len;
	}

	void gtf_assert(bool expr, const std::string& line, const std::string& msg) {
		if (!expr) {
			fprintf(stderr, "The GTF file might be corrupted!\n");
			fprintf(stderr, "Stop at line : %s\n", line.c_str());
			fprintf(stderr, "Error Message: %s\n", msg.c_str());
			exit(-1);
		}
	}  
};

#endif
