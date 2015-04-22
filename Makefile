CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-parse-alignments rsem-build-read-index rsem-run-em rsem-tbam2gbam rsem-run-gibbs rsem-calculate-credibility-intervals rsem-simulate-reads rsem-bam2wig rsem-get-unique rsem-bam2readdepth rsem-sam-validator rsem-scan-for-paired-end-reads

.PHONY : all ebseq clean

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; ${MAKE} all

rsem-extract-reference-transcripts : extractRef.cpp
	$(CC) -Wall -O3 extractRef.cpp -o rsem-extract-reference-transcripts

rsem-synthesis-reference-transcripts : synthesisRef.cpp
	$(CC) -Wall -O3 synthesisRef.cpp -o rsem-synthesis-reference-transcripts

rsem-preref : preRef.o
	$(CC) preRef.o -o rsem-preref

preRef.o : preRef.cpp
	$(CC) $(COFLAGS) preRef.cpp

rsem-parse-alignments : parseIt.o sam/libbam.a
	$(CC) -o rsem-parse-alignments parseIt.o sam/libbam.a -lz -lpthread 

parseIt.o : parseIt.cpp
	$(CC) -Wall -O2 -c -I. parseIt.cpp

rsem-build-read-index : buildReadIndex.cpp
	$(CC) -O3 buildReadIndex.cpp -o rsem-build-read-index

rsem-run-em : EM.o sam/libbam.a
	$(CC) -o rsem-run-em EM.o sam/libbam.a -lz -lpthread

EM.o : EM.cpp
	$(CC) $(COFLAGS) EM.cpp

rsem-tbam2gbam : tbam2gbam.cpp sam/libbam.a
	$(CC) -O3 -Wall tbam2gbam.cpp sam/libbam.a -lz -lpthread -o $@

rsem-bam2wig : wiggle.o sam/libbam.a bam2wig.cpp
	$(CC) -O3 -Wall bam2wig.cpp wiggle.o sam/libbam.a -lz -lpthread -o $@

rsem-bam2readdepth : wiggle.o sam/libbam.a bam2readdepth.cpp
	$(CC) -O3 -Wall bam2readdepth.cpp wiggle.o sam/libbam.a -lz -lpthread -o $@

wiggle.o: wiggle.cpp
	$(CC) $(COFLAGS) wiggle.cpp

rsem-simulate-reads : simulation.o
	$(CC) -o rsem-simulate-reads simulation.o

simulation.o : simulation.cpp
	$(CC) $(COFLAGS) simulation.cpp

rsem-run-gibbs : Gibbs.o
	$(CC) -o rsem-run-gibbs Gibbs.o -lpthread

Gibbs.o : Gibbs.cpp
	$(CC) $(COFLAGS) Gibbs.cpp

rsem-calculate-credibility-intervals : calcCI.o
	$(CC) -o rsem-calculate-credibility-intervals calcCI.o -lpthread

calcCI.o : calcCI.cpp
	$(CC) $(COFLAGS) calcCI.cpp

rsem-get-unique : getUnique.cpp sam/libbam.a
	$(CC) -O3 -Wall getUnique.cpp sam/libbam.a -lz -lpthread -o $@

rsem-sam-validator : samValidator.cpp sam/libbam.a
	$(CC) -O3 -Wall samValidator.cpp sam/libbam.a -lz -lpthread -o $@

rsem-scan-for-paired-end-reads : scanForPairedEndReads.cpp sam/libbam.a
	$(CC) -O3 -Wall scanForPairedEndReads.cpp sam/libbam.a -lz -lpthread -o $@

ebseq :
	cd EBSeq ; ${MAKE} all

clean :
	rm -f *.o *~ $(PROGRAMS)
	cd sam ; ${MAKE} clean
	cd EBSeq ; ${MAKE} clean
