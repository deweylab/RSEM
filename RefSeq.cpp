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

#include <cassert>
#include <string>
#include <fstream>

#include "RefSeq.hpp"

RefSeq::RefSeq() {
	len = 0;
	name = seq = "";
}

RefSeq::RefSeq(const std::string& name, const std::string& seq) {
   set(name, seq);
}

bool RefSeq::read(std::ifstream& fin) {
	std::string line;

	if (!getline(fin, name)) return false;
	name = name.substr(1);
	if (!getline(fin, seq)) return false;
	len = seq.length();

	return true;
}

void RefSeq::write(std::ofstream& fout) {
	fout<< ">"<< name<< std::endl;
	fout<< seq<< std::endl;
}
