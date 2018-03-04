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

#include <cstdio>
#include <string>
#include <vector>

#include "my_assert.h"
#include "GroupInfo.hpp"

void GroupInfo::load(const char* groupF) {
	FILE *fi = fopen(groupF, "r");
	int pos;

	general_assert(fi != NULL, "Cannot open " + cstrtos(groupF) + "! It may not exist.");
	
	starts.clear();
	while(fscanf(fi, "%d", &pos) == 1) {
		starts.push_back(pos);
	}
	fclose(fi);
	
	m = starts.size() - 1;
	gids = new int[starts.back()];
	for (int i = 0; i < m; ++i) {
		for (int j = starts[i]; j < starts[i + 1]; ++j) {
			gids[j] = i;
		}
	}
}
