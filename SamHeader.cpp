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

#include<cstdio>
#include<cassert>
#include<string>
#include<set>
#include<map>
#include<fstream>
#include<sstream>

#include "my_assert.h"
#include "SamHeader.hpp"

void SamHeader::replaceSQ(const char* faiF) {
  std::ifstream fin(faiF);
  general_assert(fin.is_open(), "Cannot open " + cstrtos(faiF) + "! It may not exist.");

  std::string line;
  size_t pos;
  
  SQstr = "";
  while (getline(fin, line)) {
    pos = line.find_first_of('\t');
    assert(pos != std::string::npos && pos > 0 && pos + 1 < line.length() && line[pos + 1] != '\t');
    SQstr += "@SQ\tSN:" + line.substr(0, pos) + "\tLN:" + line.substr(pos + 1, line.find_first_of('\t', pos + 1)) + "\n";
  }
  fin.close();
}

std::map<std::string, std::string> SamHeader::parse_line(const std::string& line) {
  size_t len = line.length();
  assert(line.substr(0, 3) != "@CO" && len > 4);

  size_t fr, to, colon;
  std::string field;
  std::map<std::string, std::string> dict;

  fr = 4;
  do {
    to = line.find_first_of('\t', fr);
    field = line.substr(fr, to);
    colon = field.find_first_of(':');
    if (colon != std::string::npos) {
      dict[field.substr(0, colon)] = field.substr(colon + 1);
    }
    fr = to;
  } while (fr != std::string::npos && (++fr) < len);

  return dict;
}

void SamHeader::parse_text(const char* text) {
  std::istringstream strin(text);
  std::string line, tag;

  std::map<std::string, std::string> dict;
  std::map<std::string, std::string>::iterator dict_iter;

  
  HDstr = SQstr = RGstr = PGstr = COstr = other = "";
  pids.clear();
  
  while (getline(strin, line)) {
    if (line[0] != '@') continue;
    tag = line.substr(1, 2);
    if (tag == "HD") {
      general_assert(HDstr == "", "@HD tag can only present once!");
      HDstr = line; HDstr += "\n";
    }
    else if (tag == "SQ") {
      SQstr += line; SQstr += "\n";
    }
    else if (tag == "RG") {
      RGstr += line; RGstr += "\n";
    }
    else if (tag == "PG") {
      dict = parse_line(line);
      dict_iter = dict.find("ID");
      general_assert(dict_iter != dict.end(), "\"" + line + "\" does not contain an ID!" );

      general_assert(pids.find(dict_iter->second) == pids.end(), "Program record identifier " + dict_iter->second + " is not unique!");
      pids.insert(dict_iter->second);
      
      PGstr += line; PGstr += "\n";
    }
    else if (tag == "CO") {
      COstr += line; COstr += "\n";
    }
    else {
      other += line; line += "\n";
    }
  }
}
