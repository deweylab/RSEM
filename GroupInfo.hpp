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

#ifndef GROUPINFO_H_
#define GROUPINFO_H_

#include<cassert>
#include<vector>

class GroupInfo {
public:
	GroupInfo() { m = 0; starts.clear(); gids = NULL; }
	~GroupInfo() { if (gids != NULL) delete[] gids; }

	void load(const char* groupF);

	int getm() const { return m; }

	int gidAt(int sid) const {
		assert(sid > 0 && sid < starts.back());
		return gids[sid];
	}

	// sp : start position
	int spAt(int gid) const {
		assert(gid >= 0 && gid <= m);
		return starts[gid];
	}

private:
	int m; // m genes
	std::vector<int> starts; // genes' start positions
	int *gids; // hash
};

#endif /* GROUPINFO_H_ */
