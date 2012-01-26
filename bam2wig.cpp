#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "wiggle.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: rsem-bam2wig sorted_bam_input wig_output wiggle_name\n");
		exit(-1);
	}

	UCSCWiggleTrackWriter track_writer(argv[2], argv[3]);
	build_wiggles(argv[1], track_writer);

	return 0;
}
