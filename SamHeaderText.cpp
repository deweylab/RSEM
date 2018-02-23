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

#include<cstdio>
#include<cassert>
#include<string>
#include<set>
#include<map>
#include<fstream>
#include<sstream>

#include "htslib/sam.h"
#include "my_assert.h"
#include "SamHeaderText.hpp"

void SamHeaderText::addProgram(const std::string& pid, const std::string& version, const std::string& command) {
	std::string PGstr = "@PG\tID:" + pid + "\tPN:" + pid;
	if (version != "") PGstr += "\tVN:" + version;
	if (command != "") PGstr += "\tCL:" + command;
	if (top_pid != "") PGstr += "\tPP:" + top_pid;
		
	int s = PGs.size();
	PGs.resize(s + 1);
	for (int i = s; i > 0; --i) PGs[i] = PGs[i - 1];
	PGs[0] = PGstr;

	top_pid = pid;
}

void SamHeaderText::replaceSQ(const char* faiF) {
	std::ifstream fin(faiF);
	general_assert(fin.is_open(), "Cannot open " + cstrtos(faiF) + "! It may not exist.");

	size_t pos;
	std::string line;
	
	SQs.clear();
	while (getline(fin, line)) {
		pos = line.find_first_of('\t');
		assert(pos != std::string::npos && pos > 0 && pos + 1 < line.length() && line[pos + 1] != '\t');
		SQs.push_back("@SQ\tSN:" + line.substr(0, pos) + "\tLN:" + line.substr(pos + 1, line.find_first_of('\t', pos + 1)));
	}

	fin.close();
}

bam_hdr_t* SamHeaderText::create_header() {
	std::string textstr = "";
	for (size_t i = 0; i < HDs.size(); ++i) textstr += HDs[i] + "\n";
	for (size_t i = 0; i < SQs.size(); ++i) textstr += SQs[i] + "\n";
	for (size_t i = 0; i < RGs.size(); ++i) textstr += RGs[i] + "\n";
	for (size_t i = 0; i < PGs.size(); ++i) textstr += PGs[i] + "\n";
	for (size_t i = 0; i < COs.size(); ++i) textstr += COs[i] + "\n";

	int l_text = textstr.length();
	char* text = strdup(textstr.c_str());
	bam_hdr_t *h = sam_hdr_parse(l_text, text);	
	h->l_text = l_text;
	h->text = text;

	return h;
}

std::map<std::string, std::string> SamHeaderText::parse_line(const std::string& line) {
	size_t len = line.length();
	assert(line.substr(0, 3) != "@CO" && len > 4);

	size_t fr, to, colon;
	std::string field;
	std::map<std::string, std::string> dict;

	fr = 4;
	do {
		to = line.find_first_of('\t', fr);
		field = line.substr(fr, to - fr);
		colon = field.find_first_of(':');
		if (colon != std::string::npos) {
			dict[field.substr(0, colon)] = field.substr(colon + 1);
		}
		fr = to;
	} while (fr != std::string::npos && (++fr) < len);

	return dict;
}

void SamHeaderText::parse_text(const char* text) {
	std::istringstream strin(text);
	std::string line, tag;

	HDs.clear(); SQs.clear(); RGs.clear(); PGs.clear(); COs.clear();
	while (getline(strin, line)) {
		if (line[0] != '@') continue;
		tag = line.substr(1, 2);
		if (tag == "HD") HDs.push_back(line);
		else if (tag == "SQ") SQs.push_back(line);
		else if (tag == "RG") RGs.push_back(line);
		else if (tag == "PG") PGs.push_back(line);
		else if (tag == "CO") COs.push_back(line);
		else general_assert(false, "Cannot recognize tag " + tag + "!");		
	}

	general_assert(HDs.size() <= 1, "@HD tag can only present once!");
}

void SamHeaderText::fillSQ(const bam_hdr_t* h) {
	SQs.clear();
	for (int i = 0; i < h->n_targets; ++i)
		SQs.push_back("@SQ\tSN:" + cstrtos(h->target_name[i]) + "\t" + "LN:" + itos(h->target_len[i]));
}

void SamHeaderText::assign_top_pid() {
	top_pid = "";
	if (PGs.size() > 0) {
		std::map<std::string, std::string> dict = parse_line(PGs[0]);
		std::map<std::string, std::string>::iterator it = dict.find("ID");
		general_assert(it != dict.end(), "@PG line does not have an ID tag: " + PGs[0] + "!");
		top_pid = it->second;
	}
}
