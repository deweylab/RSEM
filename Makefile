SAMTOOLS = samtools-1.3
HTSLIB = htslib-1.3

# overridable, defaulting to local copy
BOOST = .


CXX = g++
CXXFLAGS = -std=gnu++98 -Wall -I. -I$(BOOST) -I$(SAMTOOLS)/$(HTSLIB)
CPPFLAGS =

LDFLAGS =
LDLIBS =


SAMHEADERS = $(SAMTOOLS)/$(HTSLIB)/htslib/sam.h
SAMLIBS = $(SAMTOOLS)/$(HTSLIB)/libhts.a
CONFIGURE = ./configure

OBJS1 = parseIt.o
OBJS2 = extractRef.o synthesisRef.o preRef.o buildReadIndex.o wiggle.o tbam2bam.o bam2wig.o bam2readdepth.o getUnique.o samValidator.o scanForPairedEndReads.o
OBJS3 = EM.o Gibbs.o calcCI.o simulation.o

PROGS1 = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-build-read-index rsem-simulate-reads
PROGS2 = rsem-parse-alignments rsem-run-em rsem-tbam2gbam rsem-bam2wig rsem-bam2readdepth rsem-get-unique rsem-sam-validator rsem-scan-for-paired-end-reads
PROGS3 = rsem-run-gibbs rsem-calculate-credibility-intervals

PROGRAMS = $(PROGS1) $(PROGS2) $(PROGS3)


.PHONY : all ebseq clean test

all : $(PROGRAMS)

$(SAMTOOLS)/$(HTSLIB)/libhts.a : 
	cd $(SAMTOOLS) && $(CONFIGURE) --without-curses && ${MAKE} samtools


# Compile objects
$(OBJS1) :
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O2 -c -o $@ $<

$(OBJS2) :
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O3 -c -o $@ $<

$(OBJS3) :
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O3 -ffast-math -c -o $@ $<


# Generate executables
$(PROGS1) :
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(PROGS2) :
	$(CXX) $(LDFLAGS) -pthread -o $@ $^ $(LDLIBS) -lz

$(PROGS3) :
	$(CXX) $(LDFLAGS) -pthread -o $@ $^ $(LDLIBS)


# Dependencies for executables
rsem-extract-reference-transcripts : extractRef.o
rsem-synthesis-reference-transcripts : synthesisRef.o
rsem-preref : preRef.o
rsem-build-read-index : buildReadIndex.o
rsem-simulate-reads : simulation.o

rsem-parse-alignments : parseIt.o $(SAMLIBS)
rsem-run-em : EM.o $(SAMLIBS)
rsem-tbam2gbam : tbam2gbam.o $(SAMLIBS)
rsem-bam2wig : bam2wig.o wiggle.o $(SAMLIBS)
rsem-bam2readdepth : bam2readdepth.o wiggle.o $(SAMLIBS)
rsem-get-unique : getUnique.o $(SAMLIBS)
rsem-sam-validator : samValidator.o $(SAMLIBS)
rsem-scan-for-paired-end-reads : scanForPairedEndReads.o $(SAMLIBS)

rsem-run-gibbs : Gibbs.o
rsem-calculate-credibility-intervals : calcCI.o

# Dependencies for objects
parseIt.o : parseIt.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h GroupInfo.h Transcripts.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h HitContainer.h SamParser.h

extractRef.o : extractRef.cpp utils.h my_assert.h GTFItem.h Transcript.h Transcripts.h
synthesisRef.o : synthesisRef.cpp utils.h my_assert.h Transcript.h Transcripts.h
preRef.o : preRef.cpp utils.h RefSeq.h Refs.h PolyARules.h RefSeqPolicy.h AlignerRefSeqPolicy.h
buildReadIndex.o : buildReadIndex.cpp utils.h
wiggle.o: wiggle.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h wiggle.h
tbam2bam.o : tbam2gbam.cpp $(SAMHEADERS) utils.h Transcripts.h Transcript.h BamConverter.h sam_utils.h my_assert.h bc_aux.h
bam2wig.o : bam2wig.cpp utils.h my_assert.h wiggle.h
bam2readdepth.o : bam2readdepth.cpp utils.h my_assert.h wiggle.h
getUnique.o : getUnique.cpp $(SAMHEADERS) sam_utils.h utils.h 
samValidator.o : samValidator.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h
scanForPairedEndReads.o : scanForPairedEndReads.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h 

EM.o : EM.cpp $(SAMHEADERS) utils.h my_assert.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h GroupInfo.h HitContainer.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h RefSeqPolicy.h PolyARules.h Profile.h NoiseProfile.h Transcript.h Transcripts.h HitWrapper.h BamWriter.h simul.h sam_utils.h sampling.h $(BOOST)/boost/random.hpp WriteResults.h
Gibbs.o : Gibbs.cpp utils.h my_assert.h $(BOOST)/boost/random.hpp sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h 
calcCI.o : calcCI.cpp utils.h my_assert.h $(BOOST)/boost/random.hpp sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h Buffer.h 
simulation.o : simulation.cpp utils.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h RefSeq.h GroupInfo.h Transcript.h Transcripts.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h simul.h $(BOOST)/boost/random.hpp WriteResults.h

# Dependencies for header files
Transcript.h : utils.h
Transcripts.h : utils.h my_assert.h Transcript.h
BowtieRefSeqPolicy.h : RefSeqPolicy.h
RefSeq.h : utils.h
Refs.h : utils.h RefSeq.h RefSeqPolicy.h PolyARules.h
SingleRead.h : Read.h
SingleReadQ.h : Read.h
PairedEndRead.h : Read.h SingleRead.h
PairedEndReadQ.h : Read.h SingleReadQ.h
PairedEndHit.h : SingleHit.h
HitContainer.h : GroupInfo.h
sam_utils.h : $(SAMHEADERS) Transcript.h Transcripts.h
SamParser.h : $(SAMHEADERS) sam_utils.h utils.h my_assert.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Transcripts.h
simul.h : $(BOOST)/boost/random.hpp
ReadReader.h : SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h ReadIndex.h
SingleModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h SingleHit.h ReadReader.h simul.h
SingleQModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h SingleHit.h ReadReader.h simul.h
PairedEndModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h PairedEndRead.h PairedEndHit.h ReadReader.h simul.h 
PairedEndQModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h PairedEndReadQ.h PairedEndHit.h ReadReader.h simul.h
HitWrapper.h : HitContainer.h
BamWriter.h : $(SAMHEADERS) sam_utils.h utils.h my_assert.h SingleHit.h PairedEndHit.h HitWrapper.h Transcript.h Transcripts.h
sampling.h : $(BOOST)/boost/random.hpp
WriteResults.h : utils.h my_assert.h GroupInfo.h Transcript.h Transcripts.h RefSeq.h Refs.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h
bc_aux.h : $(SAMHEADERS)
BamConverter.h : $(SAMHEADERS) sam_utils.h utils.h my_assert.h bc_aux.h Transcript.h Transcripts.h
Buffer.h : my_assert.h


# Compile EBSeq
ebseq :
	cd EBSeq && ${MAKE} all

# Clean
clean :
	rm -f *.o *~ $(PROGRAMS)
	cd $(SAMTOOLS) && ${MAKE} clean
	rm -f $(SAMLIBS)
	cd EBSeq && ${MAKE} clean
