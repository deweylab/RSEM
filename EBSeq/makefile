CC = g++
PROGRAMS = EBSeq rsem-for-ebseq-calculate-clustering-info

all : $(PROGRAMS)

EBSeq : 
	./install

rsem-for-ebseq-calculate-clustering-info : calcClusteringInfo.cpp
	$(CC) -O3 -Wall calcClusteringInfo.cpp -o $@

clean : 
	rm -rf blockmodeling $(PROGRAMS) *~ 
