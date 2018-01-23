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

#include <string>
#include <sstream>

#include "SEQstring.hpp"

const char SEQstring::decode[17] = "*AC*G***T******N";
const char SEQstring::decode_r[17] = "*TG*C***A******N";

// Internal ACGTN code
const int SEQstring::codes[16] = {-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, 4};
const int SEQstring::rcodes[16] = {-1, 3, 2, -1, 1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, 4};

// toString will reset dir
std::string SEQstring::toString(char dir) {
	setDir(dir);
	std::ostringstream strout;
	for (int i = 0; i < len; ++i) strout<< baseAt(i);
	return strout.str();
}
