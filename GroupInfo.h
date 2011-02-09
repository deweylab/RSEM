#ifndef GROUPINFO_H_
#define GROUPINFO_H_

#include<cstdio>
#include<cassert>
#include<vector>

class GroupInfo {
public:
	GroupInfo() { m = 0; starts.clear(); gids = NULL; }
	~GroupInfo() { m = 0; starts.clear(); if (gids != NULL) delete[] gids; }

	void load(const char*);

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

void GroupInfo::load(const char* groupF) {
	FILE *fi = fopen(groupF, "r");
	int pos;

	if (fi == NULL) { fprintf(stderr, "Cannot open %s! It may not exist.\n", groupF); exit(-1); }

	starts.clear();
	while(fscanf(fi, "%d", &pos) == 1) {
		starts.push_back(pos);
	}
	fclose(fi);

	m = starts.size() - 1;
	gids = new int[starts.back()];
	for (int i = 0; i < m; i++) {
		for (int j = starts[i]; j < starts[i + 1]; j++) {
			gids[j] = i;
		}
	}
}

#endif /* GROUPINFO_H_ */
