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

#ifndef SAMHEADERTEXT_H_
#define SAMHEADERTEXT_H_

#include <cstdlib>
#include <string>
#include <vector>
#include <map>

#include "htslib/sam.h"

#include "my_assert.h"

class SamHeaderText {
public:
	SamHeaderText(const bam_hdr_t* h) {
		parse_text(h->text);
		if (SQs.size() == 0) fillSQ(h);
		assign_top_pid();
	}

	// return the PG ID appeared first in the text
	const std::string& getProgramID() { return top_pid; }

	void addProgram(const std::string& pid, const std::string& version = "", const std::string& command = "");

	void addComment(const std::string& comment) {
		COs.push_back("@CO\t" + comment);
	}
	
	void replaceSQ(const char* faiF);

	bam_hdr_t* create_header();

private:
	std::vector<std::string> HDs, SQs, RGs, PGs, COs;
	std::string top_pid;

	std::map<std::string, std::string> parse_line(const std::string& line);
	void parse_text(const char* text);
	void fillSQ(const bam_hdr_t* h);
	void assign_top_pid();
};

#endif
