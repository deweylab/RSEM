CC = g++
CFLAGS = -Wall -I. -I./sam -I./boost -L./sam
COFLAGS = -O3 -ffast-math -c
PROGRAMS = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-parse-alignments rsem-build-read-index rsem-run-em rsem-tbam2gbam rsem-run-gibbs rsem-calculate-credibility-intervals rsem-simulate-reads rsem-bam2wig rsem-get-unique rsem-bam2readdepth rsem-sam-validator rsem-scan-for-paired-end-reads

.PHONY : all ebseq clean

all : sam/libbam.a $(PROGRAMS)

sam/libbam.a :
	cd sam ; ${MAKE} all

ebseq :
	cd EBSeq ; ${MAKE} all


calcCI.o : calcCI.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<

EM.o : EM.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<

Gibbs.o : Gibbs.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<

preRef.o : preRef.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<

parseIt.o : parseIt.cpp
	$(CC) $(CFLAGS) -O2 -c $<

simulation.o : simulation.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<

wiggle.o: wiggle.cpp
	$(CC) $(CFLAGS) $(COFLAGS) $<


rsem-extract-reference-transcripts : extractRef.cpp
	$(CC) $(CFLAGS) -O3 $< -o $@

rsem-synthesis-reference-transcripts : synthesisRef.cpp
	$(CC) $(CFLAGS) -O3 $< -o $@

rsem-preref : preRef.o
	$(CC) $(CFLAGS) $< -o $@

rsem-parse-alignments : parseIt.o
	$(CC) $(CFLAGS) -o $@ $< -lbam -lz -lpthread

rsem-build-read-index : buildReadIndex.cpp
	$(CC) $(CFLAGS) -O3 $< -o $@

rsem-run-em : EM.o
	$(CC) $(CFLAGS) -o $@ $< -lbam -lz -lpthread

rsem-tbam2gbam : tbam2gbam.cpp
	$(CC) $(CFLAGS) -O3 $< -lbam -lz -lpthread -o $@

rsem-bam2wig : wiggle.o bam2wig.cpp
	$(CC) $(CFLAGS) -O3 $^ -lbam -lz -lpthread -o $@

rsem-bam2readdepth : wiggle.o bam2readdepth.cpp
	$(CC) $(CFLAGS) -O3 $^ -lbam -lz -lpthread -o $@

rsem-simulate-reads : simulation.o
	$(CC) $(CFLAGS) -o $@ $<

rsem-run-gibbs : Gibbs.o
	$(CC) $(CFLAGS) -o $@ $< -lpthread

rsem-calculate-credibility-intervals : calcCI.o
	$(CC) $(CFLAGS) -o $@ $< -lpthread

rsem-get-unique : getUnique.cpp
	$(CC) $(CFLAGS) -O3 $< -lbam -lz -lpthread -o $@

rsem-sam-validator : samValidator.cpp
	$(CC) $(CFLAGS) -O3 $< -lbam -lz -lpthread -o $@

rsem-scan-for-paired-end-reads : scanForPairedEndReads.cpp
	$(CC) $(CFLAGS) -O3 $< -lbam -lz -lpthread -o $@

clean :
	rm -f *.o *~ $(PROGRAMS)
	cd sam ; ${MAKE} clean
	cd EBSeq ; ${MAKE} clean
