CXX = g++
PROGRAMS = EBSeq rsem-for-ebseq-calculate-clustering-info
DEPENDENCIES = blockmodeling gplots gtools gdata caTools bitops KernSmooth

.PHONY : all EBSeq clean

all : $(PROGRAMS)

EBSeq : 
	./install

rsem-for-ebseq-calculate-clustering-info : calcClusteringInfo.cpp
	$(CXX) -O3 -Wall calcClusteringInfo.cpp -o $@

clean : 
	rm -rf $(PROGRAMS) $(DEPENDENCIES) BiocInstaller *~
