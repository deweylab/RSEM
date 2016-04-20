/* Copyright (c) 2016
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

#ifndef SAMHEADER_H_
#define SAMHEADER_H_

#include<cstdlib>
#include<string>
#include<set>
#include<map>

#include "htslib/sam.h"

class SamHeader {
public:
  SamHeader(const char* text = NULL) {
    if (text != NULL) parse_text(text);
  }

  void replaceSQ(const char* faiF);

  void insertPG(const std::string& pid, const std::string& command = "") {
    if (pids.find(pid) == pids.end()) {
      pids.insert(pid);
      PGstr += "@PG\tID:" + pid;
      if (command != "") PGstr += "\tCL:" + command;
      PGstr += "\n";
    }
  }

  void addComment(const std::string& comment) {
    COstr += "@CO\t" + comment + "\n";
  }
  
  bam_hdr_t* create_header() {
    std::string text = HDstr + SQstr + RGstr + PGstr + COstr + other;
    bam_hdr_t *h = sam_hdr_parse(text.length(), text.c_str());

    h->l_text = text.length();
    h->text = (char*)calloc(h->l_text + 1, 1);
    strcpy(h->text, text.c_str());

    return h;
  }

private:
  std::string HDstr, SQstr, RGstr, PGstr, COstr, other;
  std::set<std::string> pids;
  
  std::map<std::string, std::string> parse_line(const std::string& line);
  void parse_text(const char* text);
};

#endif
