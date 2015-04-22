CC = g++
CFLAGS = -Wall -I. -I./sam -I./boost
COFLAGS = -O3 -ffast-math -c
PROGRAMS = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-parse-alignments rsem-build-read-index rsem-run-em rsem-tbam2gbam rsem-run-gibbs rsem-calculate-credibility-intervals rsem-simulate-reads rsem-bam2wig rsem-get-unique rsem-bam2readdepth rsem-sam-validator rsem-scan-for-paired-end-reads

.PHONY : all ebseq clean

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; ${MAKE} all

rsem-extract-reference-transcripts : extractRef.cpp
	$(CC) $(CFLAGS) -O3 extractRef.cpp -o rsem-extract-reference-transcripts

rsem-synthesis-reference-transcripts : synthesisRef.cpp
	$(CC) $(CFLAGS) -O3 synthesisRef.cpp -o rsem-synthesis-reference-transcripts

rsem-preref : preRef.o
	$(CC) $(CFLAGS) preRef.o -o rsem-preref

preRef.o : preRef.cpp
	$(CC) $(CFLAGS) $(COFLAGS) preRef.cpp

rsem-parse-alignments : parseIt.o sam/libbam.a
	$(CC) $(CFLAGS) -o rsem-parse-alignments parseIt.o sam/libbam.a -lz -lpthread

parseIt.o : parseIt.cpp
	$(CC) $(CFLAGS) -O2 -c parseIt.cpp

rsem-build-read-index : buildReadIndex.cpp
	$(CC) $(CFLAGS) -O3 buildReadIndex.cpp -o rsem-build-read-index

rsem-run-em : EM.o sam/libbam.a
	$(CC) $(CFLAGS) -o rsem-run-em EM.o sam/libbam.a -lz -lpthread

EM.o : EM.cpp
	$(CC) $(CFLAGS) $(COFLAGS) EM.cpp

rsem-tbam2gbam : tbam2gbam.cpp sam/libbam.a
	$(CC) $(CFLAGS) -O3 tbam2gbam.cpp sam/libbam.a -lz -lpthread -o $@

rsem-bam2wig : wiggle.o sam/libbam.a bam2wig.cpp
	$(CC) $(CFLAGS) -O3 bam2wig.cpp wiggle.o sam/libbam.a -lz -lpthread -o $@

rsem-bam2readdepth : wiggle.o sam/libbam.a bam2readdepth.cpp
	$(CC) $(CFLAGS) -O3 bam2readdepth.cpp wiggle.o sam/libbam.a -lz -lpthread -o $@

wiggle.o: wiggle.cpp
	$(CC) $(CFLAGS) $(COFLAGS) wiggle.cpp

rsem-simulate-reads : simulation.o
	$(CC) $(CFLAGS) -o rsem-simulate-reads simulation.o

simulation.o : simulation.cpp
	$(CC) $(CFLAGS) $(COFLAGS) simulation.cpp

rsem-run-gibbs : Gibbs.o
	$(CC) $(CFLAGS) -o rsem-run-gibbs Gibbs.o -lpthread

Gibbs.o : Gibbs.cpp
	$(CC) $(CFLAGS) $(COFLAGS) Gibbs.cpp

rsem-calculate-credibility-intervals : calcCI.o
	$(CC) $(CFLAGS) -o rsem-calculate-credibility-intervals calcCI.o -lpthread

calcCI.o : calcCI.cpp
	$(CC) $(CFLAGS) $(COFLAGS) calcCI.cpp

rsem-get-unique : getUnique.cpp sam/libbam.a
	$(CC) $(CFLAGS) -O3 getUnique.cpp sam/libbam.a -lz -lpthread -o $@

rsem-sam-validator : samValidator.cpp sam/libbam.a
	$(CC) $(CFLAGS) -O3 samValidator.cpp sam/libbam.a -lz -lpthread -o $@

rsem-scan-for-paired-end-reads : scanForPairedEndReads.cpp sam/libbam.a
	$(CC) $(CFLAGS) -O3 scanForPairedEndReads.cpp sam/libbam.a -lz -lpthread -o $@

ebseq :
	cd EBSeq ; ${MAKE} all

clean :
	rm -f *.o *~ $(PROGRAMS)
	cd sam ; ${MAKE} clean
	cd EBSeq ; ${MAKE} clean
