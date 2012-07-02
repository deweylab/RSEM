#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>

#include "my_assert.h"
#include "wiggle.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: rsem-bam2readdepth sorted_bam_input readdepth_output\n");
    exit(-1);
  }

  ofstream fout(argv[2]);
  general_assert(fout.is_open(), "Cannot write to " + cstrtos(argv[2]) + "!");

  ReadDepthWriter depth_writer(fout);
  
  build_wiggles(argv[1], depth_writer);

  fout.close();

  return 0;
}
