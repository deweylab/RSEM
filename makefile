CC = g++
#LFLAGS = -Wall -O3 -ffast-math
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = rsem-bam2wig rsem-build-read-index rsem-run-em rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-parse-alignments rsem-preref rsem-simulate-reads rsem-run-gibbs rsem-calculate-credibility-intervals


all : build-sam $(PROGRAMS)

build-sam : 
	cd sam ; ${MAKE} all

Transcript.h : utils.h

Transcripts.h : Transcript.h

rsem-extract-reference-transcripts : utils.h GTFItem.h Transcript.h Transcripts.h extractRef.cpp
	$(CC) -Wall -O3 extractRef.cpp -o rsem-extract-reference-transcripts

rsem-synthesis-reference-transcripts : utils.h Transcript.h Transcripts.h synthesisRef.cpp
	$(CC) -Wall -O3 synthesisRef.cpp -o rsem-synthesis-reference-transcripts

BowtieRefSeqPolicy.h : RefSeqPolicy.h

RefSeq.h : utils.h

Refs.h : utils.h RefSeq.h RefSeqPolicy.h PolyARules.h


rsem-preref : preRef.o
	$(CC) preRef.o -o rsem-preref

preRef.o : utils.h RefSeq.h Refs.h PolyARules.h RefSeqPolicy.h AlignerRefSeqPolicy.h preRef.cpp
	$(CC) $(COFLAGS) preRef.cpp


SingleRead.h : Read.h

SingleReadQ.h : Read.h

PairedEndRead.h : Read.h SingleRead.h

PairedEndReadQ.h : Read.h SingleReadQ.h


PairedEndHit.h : SingleHit.h

HitContainer.h : GroupInfo.h


SamParser.h : sam/sam.h sam/bam.h utils.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Transcript.h Transcripts.h


rsem-parse-alignments : parseIt.o sam/libbam.a
	$(CC) -o rsem-parse-alignments parseIt.o sam/libbam.a -lz 

parseIt.o : utils.h GroupInfo.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h HitContainer.h SamParser.h Transcript.h Transcripts.h sam/sam.h sam/bam.h parseIt.cpp
	$(CC) $(COFLAGS) parseIt.cpp


rsem-build-read-index : utils.h buildReadIndex.cpp
	$(CC) -O3 buildReadIndex.cpp -o rsem-build-read-index


ReadReader.h : SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h ReadIndex.h

SingleModel.h : utils.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h SingleHit.h ReadReader.h simul.h

SingleQModel.h : utils.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h SingleHit.h ReadReader.h simul.h

PairedEndModel.h : utils.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h PairedEndRead.h PairedEndHit.h ReadReader.h simul.h 

PairedEndQModel.h : utils.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h PairedEndReadQ.h PairedEndHit.h ReadReader.h simul.h

HitWrapper.h : HitContainer.h

BamWriter.h : sam/sam.h sam/bam.h utils.h SingleHit.h PairedEndHit.h HitWrapper.h Transcript.h Transcripts.h

rsem-run-em : EM.o sam/libbam.a
	$(CC) -o rsem-run-em EM.o sam/libbam.a -lz -lpthread

EM.o : utils.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h GroupInfo.h HitContainer.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h RefSeqPolicy.h PolyARules.h Profile.h NoiseProfile.h Transcript.h Transcripts.h HitWrapper.h BamWriter.h sam/bam.h sam/sam.h EM.cpp simul.h
	$(CC) $(COFLAGS) EM.cpp

rsem-bam2wig : sam/bam.h sam/sam.h sam/libbam.a bam2wig.cpp
	$(CC) -O3 -Wall bam2wig.cpp sam/libbam.a -lz -o rsem-bam2wig

rsem-simulate-reads : simulation.o
	$(CC) -o rsem-simulate-reads simulation.o

simulation.o : utils.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h RefSeq.h GroupInfo.h Transcript.h Transcripts.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h simul.h boost/random.hpp simulation.cpp
	$(CC) $(COFLAGS) simulation.cpp

rsem-run-gibbs : Gibbs.o
	$(CC) -o rsem-run-gibbs Gibbs.o 

#some header files are omitted
Gibbs.o : utils.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h boost/random.hpp Gibbs.cpp 
	$(CC) $(COFLAGS) Gibbs.cpp

rsem-calculate-credibility-intervals : calcCI.o
	$(CC) -o rsem-calculate-credibility-intervals calcCI.o

#some header files are omitted
calcCI.o : utils.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h calcCI.cpp
	$(CC) $(COFLAGS) calcCI.cpp

clean:
	rm -f *.o *~ $(PROGRAMS)
	cd sam ; ${MAKE} clean

