#include <iostream>
#include "wiggle.h"

int main(int argc, char* argv[]) {
	if (argc != 2) {
		printf("Usage: rsem-bam2readdepth sorted_bam_input\n");
        std::exit(1);
	}
    ReadDepthWriter depth_writer(std::cout);
    build_wiggles(argv[1], depth_writer);

	return 0;
}
