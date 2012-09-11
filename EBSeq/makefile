CC = g++
PROGRAMS = blockmodeling EBSeq rsem-for-ebseq-calculate-clustering-info

all : $(PROGRAMS)

blockmodeling : blockmodeling_0.1.8.tar.gz
	R CMD INSTALL -l "." blockmodeling_0.1.8.tar.gz

EBSeq : blockmodeling EBSeq_1.1.3.tar.gz
	R CMD INSTALL -l "." EBSeq_1.1.3.tar.gz

rsem-for-ebseq-calculate-clustering-info : calcClusteringInfo.cpp
	$(CC) -O3 -Wall calcClusteringInfo.cpp -o $@

clean : 
	rm -rf $(PROGRAMS) *~
