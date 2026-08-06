SHELL := /bin/bash

SAMTOOLS = samtools-1.3
HTSLIB = htslib-1.3

ifneq ($(cygwin), true)
  SAMTOOLS_MAKEFILE = Makefile
else
  SAMTOOLS_MAKEFILE = Makefile.cygwin
endif

# Compilation variables
CXX = g++
CXXFLAGS = -std=c++17 -Wall -I. -I$(SAMTOOLS)/$(HTSLIB)
CPPFLAGS =

LDFLAGS =
LDLIBS =

# Installation variables
INSTALL = install
INSTALL_PROGRAM = $(INSTALL) -p
INSTALL_DATA = $(INSTALL) -p -m 644
INSTALL_DIR = $(INSTALL) -d
STRIP ?=strip

prefix ?= /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin

# Auxiliary variables for compilation
SAMHEADERS = $(SAMTOOLS)/$(HTSLIB)/htslib/sam.h
SAMLIBS = $(SAMTOOLS)/$(HTSLIB)/libhts.a
CONFIGURE = ./configure

OBJS1 = parseIt.o
OBJS2 = extractRef.o synthesisRef.o preRef.o buildReadIndex.o wiggle.o tbam2gbam.o bam2wig.o bam2readdepth.o getUnique.o samValidator.o scanForPairedEndReads.o SamHeader.o
OBJS3 = EM.o Gibbs.o calcCI.o simulation.o

PROGS1 = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-build-read-index rsem-simulate-reads
PROGS2 = rsem-parse-alignments rsem-run-em rsem-tbam2gbam rsem-bam2wig rsem-bam2readdepth rsem-get-unique rsem-sam-validator rsem-scan-for-paired-end-reads
PROGS3 = rsem-run-gibbs rsem-calculate-credibility-intervals

PROGRAMS = $(PROGS1) $(PROGS2) $(PROGS3)

# Auxiliary variables for installation
SCRIPTS = rsem-prepare-reference rsem-calculate-expression rsem-refseq-extract-primary-assembly rsem-gff3-to-gtf rsem-plot-model \
	  rsem-plot-transcript-wiggles rsem-gen-transcript-plots rsem-generate-data-matrix \
	  extract-transcript-to-gene-map-from-trinity convert-sam-for-rsem    




.PHONY : all ebseq pRSEM clean

all : $(PROGRAMS) $(SAMTOOLS)/samtools

$(SAMTOOLS)/samtools :
	cd $(SAMTOOLS) && $(CONFIGURE) --without-curses && $(MAKE) -f $(SAMTOOLS_MAKEFILE) samtools

$(SAMLIBS) : $(SAMTOOLS)/samtools


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
rsem-run-em : EM.o SamHeader.o $(SAMLIBS)
rsem-tbam2gbam : tbam2gbam.o SamHeader.o $(SAMLIBS)
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
tbam2gbam.o : tbam2gbam.cpp $(SAMHEADERS) utils.h Transcripts.h Transcript.h BamConverter.h sam_utils.h SamHeader.hpp my_assert.h bc_aux.h
bam2wig.o : bam2wig.cpp utils.h my_assert.h wiggle.h
bam2readdepth.o : bam2readdepth.cpp utils.h my_assert.h wiggle.h
getUnique.o : getUnique.cpp $(SAMHEADERS) sam_utils.h utils.h 
samValidator.o : samValidator.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h
scanForPairedEndReads.o : scanForPairedEndReads.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h 
SamHeader.o : SamHeader.cpp $(SAMHEADERS) SamHeader.hpp 

EM.o : EM.cpp $(SAMHEADERS) utils.h my_assert.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h GroupInfo.h HitContainer.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h RefSeqPolicy.h PolyARules.h Profile.h NoiseProfile.h Transcript.h Transcripts.h HitWrapper.h BamWriter.h simul.h sam_utils.h SamHeader.hpp sampling.h WriteResults.h
Gibbs.o : Gibbs.cpp utils.h my_assert.h sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h 
calcCI.o : calcCI.cpp utils.h my_assert.h sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h Buffer.h 
simulation.o : simulation.cpp utils.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h RefSeq.h GroupInfo.h Transcript.h Transcripts.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h simul.h WriteResults.h

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
simul.h :
ReadReader.h : SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h ReadIndex.h
SingleModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h SingleHit.h ReadReader.h simul.h
SingleQModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h SingleHit.h ReadReader.h simul.h
PairedEndModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h Profile.h NoiseProfile.h ModelParams.h RefSeq.h Refs.h SingleRead.h PairedEndRead.h PairedEndHit.h ReadReader.h simul.h 
PairedEndQModel.h : utils.h my_assert.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h Refs.h SingleReadQ.h PairedEndReadQ.h PairedEndHit.h ReadReader.h simul.h
HitWrapper.h : HitContainer.h
BamWriter.h : $(SAMHEADERS) sam_utils.h SamHeader.hpp utils.h my_assert.h SingleHit.h PairedEndHit.h HitWrapper.h Transcript.h Transcripts.h
sampling.h :
WriteResults.h : utils.h my_assert.h GroupInfo.h Transcript.h Transcripts.h RefSeq.h Refs.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h
bc_aux.h : $(SAMHEADERS)
BamConverter.h : $(SAMHEADERS) sam_utils.h SamHeader.hpp utils.h my_assert.h bc_aux.h Transcript.h Transcripts.h
Buffer.h : my_assert.h
SamHeader.hpp : $(SAMHEADERS)

# Compile EBSeq
ebseq :
	cd EBSeq && $(MAKE) all

# Compile pRSEM
pRSEM :
	cd pRSEM && $(MAKE) all


# Install RSEM
install : $(PROGRAMS) $(SCRIPTS) $(SAMTOOLS)/samtools rsem_utils.py
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(bindir)/$(SAMTOOLS)
	$(foreach prog,$(PROGRAMS),$(INSTALL_PROGRAM) $(prog) $(DESTDIR)$(bindir)/$(prog) ; $(STRIP) $(DESTDIR)$(bindir)/$(prog) ;)
	$(INSTALL_PROGRAM) $(SAMTOOLS)/samtools $(DESTDIR)$(bindir)/$(SAMTOOLS)/samtools
	$(STRIP) $(DESTDIR)$(bindir)/$(SAMTOOLS)/samtools
	$(foreach script,$(SCRIPTS),$(INSTALL_PROGRAM) $(script) $(DESTDIR)$(bindir)/$(script) ;)
	$(INSTALL_DATA) rsem_utils.py $(DESTDIR)$(bindir)/rsem_utils.py

# Clean
clean :
	rm -f *.o *~ $(PROGRAMS)
	cd $(SAMTOOLS) && $(MAKE) clean-all
	cd EBSeq && $(MAKE) clean
	cd pRSEM && $(MAKE) clean

# ===================================================================
# ---- Test targets -------------------------------------------------
# ===================================================================

# Gold layout: tests/gold/<read_mode>/<aligner>/{reference,expression,expression_ci,simulated}
# read_mode: single_end (reads_1.fastq) | paired_end (reads_1 + reads_2, --paired-end)
READ_MODES := single_end paired_end
ALIGNERS := bowtie bowtie2 hisat2 star

# STAR 2.7.6a for gold + star tests (Bioconda binary + wrapper; run `make fetch-star-276a` first).
STAR_276A_DIR := $(CURDIR)/tests/tools/star-2.7.6a
STAR_EXTRA_ARGS_bowtie :=
STAR_EXTRA_ARGS_bowtie2 :=
STAR_EXTRA_ARGS_hisat2 :=
STAR_EXTRA_ARGS_star := --star-path $(STAR_276A_DIR)

# STAR --genomeSAindexNbases for the small viral test reference (stable SAindex size + reproducible gold).
STAR_GENOME_SA_INDEX_NBASES := 6

# Exclude machine-specific files from the reference diff:
#   Log.out            - STAR timestamp
#   genomeParameters.txt - embeds STAR's absolute path; diffed separately below
#   *.ht2              - HISAT2 binary index, not byte-reproducible (covered by expression tests)
REF_DIFF_EXCL_bowtie :=
REF_DIFF_EXCL_bowtie2 :=
REF_DIFF_EXCL_hisat2 := -x '*.ht2'
REF_DIFF_EXCL_star := -x Log.out -x genomeParameters.txt

TEST_TARGETS := $(foreach m,$(READ_MODES),$(foreach a,$(ALIGNERS),test-prepare-reference-$(m)-$(a) test-calculate-expression-$(m)-$(a) test-calculate-expression-ci-$(m)-$(a) test-simulate-reads-$(m)-$(a)))

.PHONY: test test-all generate-gold fetch-star-276a $(TEST_TARGETS) \
	test-prepare-reference test-calculate-expression test-calculate-expression-ci test-simulate-reads

# Convenience aliases (all read modes × aligners)
test-prepare-reference: $(foreach m,$(READ_MODES),$(foreach a,$(ALIGNERS),test-prepare-reference-$(m)-$(a)))
test-calculate-expression: $(foreach m,$(READ_MODES),$(foreach a,$(ALIGNERS),test-calculate-expression-$(m)-$(a)))
test-calculate-expression-ci: $(foreach m,$(READ_MODES),$(foreach a,$(ALIGNERS),test-calculate-expression-ci-$(m)-$(a)))
test-simulate-reads: $(foreach m,$(READ_MODES),$(foreach a,$(ALIGNERS),test-simulate-reads-$(m)-$(a)))

# Default: fail fast. Depends on 'all' to ensure binaries are built first.
# fetch-star-276a is a no-op if STAR 2.7.6a is already unpacked (required for star/* tests).
test: all fetch-star-276a $(TEST_TARGETS) test-options

# CI / full signal: keep going
test-all: all
	$(MAKE) -k test

# Print every individual test target (matrix + option cases), space-separated.
# Used by CI to drive each target separately and build a per-target result table,
# without duplicating the matrix definition in the workflow.
.PHONY: print-test-targets
print-test-targets:
	@echo $(TEST_TARGETS) $(foreach _c,$(OPTION_CASES),test-option-$(_c))

# Per-aligner flags for rsem-prepare-reference / rsem-calculate-expression
PREP_FLAGS_bowtie := --bowtie
PREP_FLAGS_bowtie2 := --bowtie2
PREP_FLAGS_hisat2 := --hisat2-hca
# Uses STAR_GENOME_SA_INDEX_NBASES (see rsem-prepare-reference --star-genome-sa-index-nbases).
PREP_FLAGS_star := --star --star-genome-sa-index-nbases $(STAR_GENOME_SA_INDEX_NBASES)
# Bowtie 1 is the default aligner when no --bowtie2/--star/--hisat2-hca is passed
CALC_FLAGS_bowtie :=
CALC_FLAGS_bowtie2 := --bowtie2
CALC_FLAGS_hisat2 := --hisat2-hca
CALC_FLAGS_star := --star

# Generate gold standard: run pipeline per read_mode × aligner under tests/gold/<read_mode>/<aligner>/
# Requires bowtie, bowtie2, hisat2, and STAR 2.7.6a under tests/tools/star-2.7.6a/ (see fetch-star-276a).
generate-gold: all fetch-star-276a
	@echo "==> Generating gold standard (read_modes=$(READ_MODES), aligners=$(ALIGNERS), TEST_THREADS=$(TEST_THREADS), TEST_SEED=$(TEST_SEED))"
	rm -rf $(TEST_OUTPUT) $(GOLD_ROOT)
	@set -e; for m in $(READ_MODES); do \
	  for a in $(ALIGNERS); do \
	    echo "==> Gold read_mode=$$m aligner=$$a"; \
	    G="$(GOLD_ROOT)/$$m/$$a"; \
	    O="$(TEST_OUTPUT)/$$m/$$a"; \
	    mkdir -p "$$O/reference" "$$O/expression/$(SAMPLE_NAME).stat" "$$O/simulated" "$$O/expression_ci/$(SAMPLE_NAME_CI).stat" \
	      "$$G/reference" "$$G/expression/$(SAMPLE_NAME).stat" "$$G/simulated" "$$G/expression_ci/$(SAMPLE_NAME_CI).stat"; \
	    starp=""; test "$$a" != star || starp="--star-path $(STAR_276A_DIR)"; \
	    case $$a in \
	      bowtie) prepflags="--bowtie"; calcflags="" ;; \
	      bowtie2) prepflags="--bowtie2"; calcflags="--bowtie2" ;; \
	      hisat2) prepflags="--hisat2-hca"; calcflags="--hisat2-hca" ;; \
	      star) prepflags="--star --star-genome-sa-index-nbases $(STAR_GENOME_SA_INDEX_NBASES)"; calcflags="--star" ;; \
	      *) echo "Unknown aligner $$a"; exit 1 ;; \
	    esac; \
	    ./rsem-prepare-reference --gtf $(TEST_GTF) $$prepflags $$starp -p $(TEST_THREADS) $(TEST_GENOME) "$$O/reference/$(REF_NAME)"; \
	    cp -a "$$O/reference/." "$$G/reference/"; \
	    if [ "$$m" = paired_end ]; then \
	      ./rsem-calculate-expression $$calcflags $$starp -p $(TEST_THREADS) --seed $(TEST_SEED) --paired-end $(TEST_READS_1) $(TEST_READS_2) "$$O/reference/$(REF_NAME)" "$$O/expression/$(SAMPLE_NAME)"; \
	    else \
	      ./rsem-calculate-expression $$calcflags $$starp -p $(TEST_THREADS) --seed $(TEST_SEED) $(TEST_READS_1) "$$O/reference/$(REF_NAME)" "$$O/expression/$(SAMPLE_NAME)"; \
	    fi; \
	    cp "$$O/expression/$(SAMPLE_NAME).genes.results" "$$O/expression/$(SAMPLE_NAME).isoforms.results" "$$G/expression/"; \
	    cp "$$O/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt" "$$O/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model" "$$O/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta" "$$G/expression/$(SAMPLE_NAME).stat/"; \
	    theta0=$$(awk 'NR==3 {print $$1}' "$$O/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta"); \
	    ./rsem-simulate-reads "$$O/reference/$(REF_NAME)" "$$O/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model" "$$O/expression/$(SAMPLE_NAME).isoforms.results" "$$theta0" 1000 "$$O/simulated/$(SAMPLE_NAME).simulated" --seed $(TEST_SEED); \
	    if [ "$$m" = paired_end ]; then \
	      cp "$$O/simulated/$(SAMPLE_NAME).simulated_1.fq" "$$O/simulated/$(SAMPLE_NAME).simulated_2.fq" "$$O/simulated/$(SAMPLE_NAME).simulated.sim.isoforms.results" "$$O/simulated/$(SAMPLE_NAME).simulated.sim.genes.results" "$$G/simulated/"; \
	    else \
	      cp "$$O/simulated/$(SAMPLE_NAME).simulated.fq" "$$O/simulated/$(SAMPLE_NAME).simulated.sim.isoforms.results" "$$O/simulated/$(SAMPLE_NAME).simulated.sim.genes.results" "$$G/simulated/"; \
	    fi; \
	    if [ "$$m" = paired_end ]; then \
	      ./rsem-calculate-expression $$calcflags $$starp -p $(TEST_THREADS) --seed $(TEST_SEED) --calc-ci --paired-end $(TEST_READS_1) $(TEST_READS_2) "$$G/reference/$(REF_NAME)" "$$O/expression_ci/$(SAMPLE_NAME_CI)"; \
	    else \
	      ./rsem-calculate-expression $$calcflags $$starp -p $(TEST_THREADS) --seed $(TEST_SEED) --calc-ci $(TEST_READS_1) "$$G/reference/$(REF_NAME)" "$$O/expression_ci/$(SAMPLE_NAME_CI)"; \
	    fi; \
	    cp "$$O/expression_ci/$(SAMPLE_NAME_CI).genes.results" "$$O/expression_ci/$(SAMPLE_NAME_CI).isoforms.results" "$$G/expression_ci/"; \
	    cp "$$O/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).cnt" "$$O/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).model" "$$O/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).theta" "$$G/expression_ci/$(SAMPLE_NAME_CI).stat/"; \
	  done; \
	done
	@echo "==> Gold standard generated under $(GOLD_ROOT)/{$(READ_MODES)}/{$(ALIGNERS)}/"

fetch-star-276a:
	@chmod +x tests/fetch-star-276a.sh 2>/dev/null || true
	tests/fetch-star-276a.sh
	@chmod +x $(STAR_276A_DIR)/STAR

# ---- Paths --------------------------------------------------------

TEST_DATA    := tests/data
TEST_OUTPUT  := tests/output
GOLD_ROOT    := tests/gold
REF_NAME     := my_ref
SAMPLE_NAME  := my_sample
SAMPLE_NAME_CI := my_sample_ci

# Reproducibility: same thread count and seed for all test commands (single thread for deterministic Gibbs/calcCI)
TEST_THREADS := 1
TEST_SEED    := 0

# Input files (SARS-CoV-2 reference + SRR11550043 read subset)
TEST_GENOME  := $(TEST_DATA)/sarscov2.fasta
TEST_GTF     := $(TEST_DATA)/sarscov2.gtf
TEST_READS_1 := $(TEST_DATA)/reads_1.fastq
TEST_READS_2 := $(TEST_DATA)/reads_2.fastq

# Per-read-mode arguments for rsem-calculate-expression
CALC_READ_ARGS_single_end := $(TEST_READS_1)
CALC_READ_ARGS_paired_end := --paired-end $(TEST_READS_1) $(TEST_READS_2)

# ---- Individual tests (read_mode × aligner) ---

define _RULE_TEST_PREP_REF
test-prepare-reference-$(2)-$(1): $$(if $$(filter star,$(1)),fetch-star-276a)
	@echo "==> Testing rsem-prepare-reference ($(2), $(1))"
	rm -rf $(TEST_OUTPUT)/$(2)/$(1)/reference
	mkdir -p $(TEST_OUTPUT)/$(2)/$(1)/reference
	./rsem-prepare-reference \
		--gtf $(TEST_GTF) \
		$(PREP_FLAGS_$(1)) \
		$(STAR_EXTRA_ARGS_$(1)) \
		-p $(TEST_THREADS) \
		$(TEST_GENOME) \
		$(TEST_OUTPUT)/$(2)/$(1)/reference/$(REF_NAME)
	diff -r $(REF_DIFF_EXCL_$(1)) $(TEST_OUTPUT)/$(2)/$(1)/reference $(GOLD_ROOT)/$(2)/$(1)/reference
	$$(if $$(filter star,$(1)),diff <(grep -v '^###' $(TEST_OUTPUT)/$(2)/$(1)/reference/genomeParameters.txt) <(grep -v '^###' $(GOLD_ROOT)/$(2)/$(1)/reference/genomeParameters.txt))
	@echo "==> test-prepare-reference-$(2)-$(1): OK"
endef
$(foreach _m,$(READ_MODES),$(foreach _a,$(ALIGNERS),$(eval $(call _RULE_TEST_PREP_REF,$(_a),$(_m)))))

define _RULE_TEST_CALC_EXPR
test-calculate-expression-$(2)-$(1): $$(if $$(filter star,$(1)),fetch-star-276a)
	@echo "==> Testing rsem-calculate-expression ($(2), $(1))"
	rm -rf $(TEST_OUTPUT)/$(2)/$(1)/expression
	mkdir -p $(TEST_OUTPUT)/$(2)/$(1)/expression
	./rsem-calculate-expression \
		$(CALC_FLAGS_$(1)) \
		$(STAR_EXTRA_ARGS_$(1)) \
		-p $(TEST_THREADS) \
		--seed $(TEST_SEED) \
		$(CALC_READ_ARGS_$(2)) \
		$(GOLD_ROOT)/$(2)/$(1)/reference/$(REF_NAME) \
		$(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME)
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME).genes.results $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).genes.results
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME).isoforms.results $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).isoforms.results
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt
	python3 tests/compare_floats.py --exact $(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model
	python3 tests/compare_floats.py $(TEST_OUTPUT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta
	@echo "==> test-calculate-expression-$(2)-$(1): OK"
endef
$(foreach _m,$(READ_MODES),$(foreach _a,$(ALIGNERS),$(eval $(call _RULE_TEST_CALC_EXPR,$(_a),$(_m)))))

define _RULE_TEST_CALC_CI
test-calculate-expression-ci-$(2)-$(1): $$(if $$(filter star,$(1)),fetch-star-276a)
	@echo "==> Testing rsem-calculate-expression with --calc-ci ($(2), $(1))"
	rm -rf $(TEST_OUTPUT)/$(2)/$(1)/expression_ci
	mkdir -p $(TEST_OUTPUT)/$(2)/$(1)/expression_ci
	./rsem-calculate-expression \
		$(CALC_FLAGS_$(1)) \
		$(STAR_EXTRA_ARGS_$(1)) \
		-p $(TEST_THREADS) \
		--seed $(TEST_SEED) \
		--calc-ci \
		$(CALC_READ_ARGS_$(2)) \
		$(GOLD_ROOT)/$(2)/$(1)/reference/$(REF_NAME) \
		$(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI)
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).genes.results $(GOLD_ROOT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).genes.results
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).isoforms.results $(GOLD_ROOT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).isoforms.results
	diff $(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).cnt $(GOLD_ROOT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).cnt
	python3 tests/compare_floats.py --exact $(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).model $(GOLD_ROOT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).model
	python3 tests/compare_floats.py $(TEST_OUTPUT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).theta $(GOLD_ROOT)/$(2)/$(1)/expression_ci/$(SAMPLE_NAME_CI).stat/$(SAMPLE_NAME_CI).theta
	@echo "==> test-calculate-expression-ci-$(2)-$(1): OK"
endef
$(foreach _m,$(READ_MODES),$(foreach _a,$(ALIGNERS),$(eval $(call _RULE_TEST_CALC_CI,$(_a),$(_m)))))

define _RULE_TEST_SIM
test-simulate-reads-$(2)-$(1):
	@echo "==> Testing rsem-simulate-reads ($(2), $(1))"
	rm -rf $(TEST_OUTPUT)/$(2)/$(1)/simulated
	mkdir -p $(TEST_OUTPUT)/$(2)/$(1)/simulated
	@theta0=`awk 'NR==3 {print $$$$1}' $(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta` && \
	./rsem-simulate-reads \
		$(GOLD_ROOT)/$(2)/$(1)/reference/$(REF_NAME) \
		$(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model \
		$(GOLD_ROOT)/$(2)/$(1)/expression/$(SAMPLE_NAME).isoforms.results \
		$$$$theta0 \
		1000 \
		$(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated \
		--seed $(TEST_SEED)
	@if [ "$(2)" = "paired_end" ]; then \
	  diff $(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated_1.fq $(GOLD_ROOT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated_1.fq; \
	  diff $(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated_2.fq $(GOLD_ROOT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated_2.fq; \
	else \
	  diff $(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.fq $(GOLD_ROOT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.fq; \
	fi
	diff $(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.sim.isoforms.results $(GOLD_ROOT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.sim.isoforms.results
	diff $(TEST_OUTPUT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.sim.genes.results $(GOLD_ROOT)/$(2)/$(1)/simulated/$(SAMPLE_NAME).simulated.sim.genes.results
	@echo "==> test-simulate-reads-$(2)-$(1): OK"
endef
$(foreach _m,$(READ_MODES),$(foreach _a,$(ALIGNERS),$(eval $(call _RULE_TEST_SIM,$(_a),$(_m)))))

# ===================================================================
# ---- Option-coverage tests -----------------------------------------
# ===================================================================
# Each case exercises a cluster of rsem-calculate-expression flags.
# tests/check_option_coverage.py verifies each flag changes the output.

OPTION_CASES := bowtie_custom bowtie2_custom gibbs_sampling fragment_modeling star_gzip_genomebam hisat2_path

# Dir passed to --hisat2-path; holds a probe wrapper that records invocations
# and delegates by absolute path. Kept off PATH so only the flag can reach it.
HISAT2_PROBE_DIR := $(CURDIR)/tests/tools/hisat2-probe
HISAT2_REAL      := $(shell command -v hisat2 2>/dev/null)
HISAT2_PROBE_LOG := $(TEST_OUTPUT)/options/hisat2_path/hisat2-probe.log

TEST_READS_1_GZ := $(TEST_DATA)/reads_1.fastq.gz
TEST_READS_2_GZ := $(TEST_DATA)/reads_2.fastq.gz

$(TEST_READS_1_GZ): $(TEST_READS_1)
	gzip -k -f $<

$(TEST_READS_2_GZ): $(TEST_READS_2)
	gzip -k -f $<

OPT_ARGS_bowtie_custom       := --bowtie-n 3 --bowtie-e 200 --bowtie-m 5 --seed-length 28 --paired-end $(TEST_READS_1) $(TEST_READS_2)
OPT_ARGS_bowtie2_custom      := --bowtie2 --bowtie2-mismatch-rate 0.05 --bowtie2-k 5 --bowtie2-sensitivity-level very_fast --paired-end $(TEST_READS_1) $(TEST_READS_2)
OPT_ARGS_gibbs_sampling      := --seed $(TEST_SEED) --paired-end --single-cell-prior --calc-pme --calc-ci --gibbs-burnin 20 --gibbs-number-of-samples 800 --gibbs-sampling-gap 2 --ci-credibility-level 0.80 --ci-number-of-samples-per-count-vector 30 $(TEST_READS_1) $(TEST_READS_2)
OPT_ARGS_fragment_modeling   := --fragment-length-mean 400 --fragment-length-sd 50 --fragment-length-min 100 --fragment-length-max 700 --estimate-rspd --num-rspd-bins 40 $(TEST_READS_1)
OPT_ARGS_star_gzip_genomebam := --star --star-path $(STAR_276A_DIR) --star-gzipped-read-file --star-output-genome-bam --paired-end $(TEST_READS_1_GZ) $(TEST_READS_2_GZ)
OPT_ARGS_hisat2_path         := --hisat2-hca --hisat2-path $(HISAT2_PROBE_DIR) --seed $(TEST_SEED) --paired-end $(TEST_READS_1) $(TEST_READS_2)

# Reference for each case (bowtie is the default when no aligner flag is given).
OPT_REF_bowtie_custom       := $(GOLD_ROOT)/paired_end/bowtie/reference/$(REF_NAME)
OPT_REF_bowtie2_custom      := $(GOLD_ROOT)/paired_end/bowtie2/reference/$(REF_NAME)
OPT_REF_gibbs_sampling      := $(GOLD_ROOT)/paired_end/bowtie/reference/$(REF_NAME)
OPT_REF_fragment_modeling   := $(GOLD_ROOT)/single_end/bowtie/reference/$(REF_NAME)
OPT_REF_star_gzip_genomebam := $(GOLD_ROOT)/paired_end/star/reference/$(REF_NAME)
OPT_REF_hisat2_path         := $(GOLD_ROOT)/paired_end/hisat2/reference/$(REF_NAME)

OPT_PREREQ_bowtie_custom       :=
OPT_PREREQ_bowtie2_custom      :=
OPT_PREREQ_gibbs_sampling      :=
OPT_PREREQ_fragment_modeling   :=
OPT_PREREQ_star_gzip_genomebam := fetch-star-276a $(TEST_READS_1_GZ) $(TEST_READS_2_GZ)
OPT_PREREQ_hisat2_path         :=

# Optional per-case environment prefix, and assertion run after the gold diffs
# for checks gold files cannot express. Unset cases expand to nothing.
OPT_ENV_hisat2_path  := RSEM_TEST_HISAT2_PROBE=$(HISAT2_PROBE_LOG) RSEM_TEST_HISAT2_BIN=$(HISAT2_REAL)
OPT_POST_hisat2_path := tests/assert_path_flag_honored.sh $(HISAT2_PROBE_LOG) $(HISAT2_PROBE_DIR) hisat2 --hisat2-path

.PHONY: test-options generate-gold-options $(foreach _c,$(OPTION_CASES),test-option-$(_c))

define _RULE_TEST_OPTION
test-option-$(1): $(OPT_PREREQ_$(1))
	@echo "==> Testing option case: $(1)"
	rm -rf $(TEST_OUTPUT)/options/$(1)
	mkdir -p $(TEST_OUTPUT)/options/$(1)
	$$(OPT_ENV_$(1)) PATH="$(CURDIR)/tests/shims:$$$$PATH" ./rsem-calculate-expression -p $(TEST_THREADS) $$(OPT_ARGS_$(1)) $$(OPT_REF_$(1)) $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME)
	diff $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).genes.results $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).genes.results
	diff $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).isoforms.results $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).isoforms.results
	diff $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt
	python3 tests/compare_floats.py --exact $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model
	python3 tests/compare_floats.py $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta
	$$(OPT_POST_$(1))
	@echo "==> test-option-$(1): OK"
endef
$(foreach _c,$(OPTION_CASES),$(eval $(call _RULE_TEST_OPTION,$(_c))))

test-options: $(foreach _c,$(OPTION_CASES),test-option-$(_c))

# Generate gold for the option-coverage cases. One rule per case via $(eval): the
# OPT_ARGS_$(1)/OPT_REF_$(1) lookups need a literal $(1), so a shell loop won't work.
.PHONY: $(foreach _c,$(OPTION_CASES),generate-gold-option-$(_c))

define _RULE_GENERATE_GOLD_OPTION
generate-gold-option-$(1): all $(OPT_PREREQ_$(1))
	@echo "==> Gold option case=$(1)"
	mkdir -p $(TEST_OUTPUT)/options/$(1) $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).stat
	$(OPT_ENV_$(1)) PATH="$(CURDIR)/tests/shims:$$$$PATH" ./rsem-calculate-expression -p $(TEST_THREADS) $(OPT_ARGS_$(1)) $(OPT_REF_$(1)) $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME)
	cp $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).genes.results $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).isoforms.results $(GOLD_ROOT)/options/$(1)/
	cp $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).cnt $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).model $(TEST_OUTPUT)/options/$(1)/$(SAMPLE_NAME).stat/$(SAMPLE_NAME).theta $(GOLD_ROOT)/options/$(1)/$(SAMPLE_NAME).stat/
endef
$(foreach _c,$(OPTION_CASES),$(eval $(call _RULE_GENERATE_GOLD_OPTION,$(_c))))

generate-gold-options: $(foreach _c,$(OPTION_CASES),generate-gold-option-$(_c))
	@echo "==> Gold option-coverage standard generated under $(GOLD_ROOT)/options/{$(OPTION_CASES)}/"
