#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<iostream>
#include<fstream>
#include<vector>

#include "utils.h"
#include "my_assert.h"
#include "Read.h"
#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Refs.h"
#include "GroupInfo.h"
#include "Transcript.h"
#include "Transcripts.h"

#include "simul.h"

using namespace std;

const int OFFSITE = 5;

READ_INT_TYPE N;
int model_type, M, m;

Refs refs;
GroupInfo gi;
Transcripts transcripts;

vector<double> eel;
vector<double> theta, counts;

int n_os;
ostream *os[2];
char outReadF[2][STRLEN];

char refF[STRLEN], groupF[STRLEN], tiF[STRLEN];
char geneResF[STRLEN], isoResF[STRLEN];

simul sampler;

void genOutReadStreams(int type, char *outFN) {
	switch(type) {
	case 0 :
		n_os = 1;
		sprintf(outReadF[0], "%s.fa", outFN);
		break;
	case 1 :
		n_os = 1;
		sprintf(outReadF[0], "%s.fq", outFN);
		break;
	case 2 :
		n_os = 2;
		for (int i = 0; i < n_os; i++)
			sprintf(outReadF[i], "%s_%d.fa", outFN, i + 1);
		break;
	case 3 :
		n_os = 2;
		for (int i = 0; i < n_os; i++)
			sprintf(outReadF[i], "%s_%d.fq", outFN, i + 1);
		break;
	}

	for (int i = 0; i < n_os; i++)
		os[i] = new ofstream(outReadF[i]);
}

template<class ModelType>
void calcExpectedEffectiveLengths(ModelType& model) {
	int lb, ub, span;
	double *pdf = NULL, *cdf = NULL, *clen = NULL; // clen[i] = sigma_{j=1}^{i}pdf[i]*(lb+i)
  
	model.getGLD().copyTo(pdf, cdf, lb, ub, span);
	clen = new double[span + 1];
	clen[0] = 0.0;
	for (int i = 1; i <= span; i++) {
		clen[i] = clen[i - 1] + pdf[i] * (lb + i);
	}

	eel.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++) {
		int totLen = refs.getRef(i).getTotLen();
		int fullLen = refs.getRef(i).getFullLen();
		int pos1 = max(min(totLen - fullLen + 1, ub) - lb, 0);
		int pos2 = max(min(totLen, ub) - lb, 0);

		if (pos2 == 0) { eel[i] = 0.0; continue; }
    
		eel[i] = fullLen * cdf[pos1] + ((cdf[pos2] - cdf[pos1]) * (totLen + 1) - (clen[pos2] - clen[pos1]));
		assert(eel[i] >= 0);
		if (eel[i] < MINEEL) { eel[i] = 0.0; }
	}
  
	delete[] pdf;
	delete[] cdf;
	delete[] clen;
}

template<class ReadType, class ModelType>
void simulate(char* modelF, char* resultsF) {
	ModelType model(&refs);
	ReadType read;
	int sid;

	model.read(modelF);
	
	//calculate eel
	calcExpectedEffectiveLengths<ModelType>(model);

	//generate theta vector
	ifstream fin(resultsF);
	string line;
	double tpm;
	double denom = 0.0;
	getline(fin, line); // read the first line, which is just column names
	for (int i = 1; i <= M; i++) {
	  getline(fin, line);
	  size_t pos = 0;
	  for (int j = 0; j < OFFSITE; j++) pos = line.find_first_of('\t', pos) + 1;
	  size_t pos2 = line.find_first_of('\t', pos);
	  if (pos2 == string::npos) pos2 = line.length();
	  tpm = atof(line.substr(pos, pos2 - pos).c_str());
	  theta[i] = tpm * eel[i]; // during simulation, there is no check for effL < 0. The reason is for that case, eel[i] here = 0 and therefore no chance to sample from it
	  denom += theta[i];
	}
	assert(denom > EPSILON);
	fin.close();
	for (int i = 1; i <= M; i++) theta[i] = theta[i] / denom * (1.0 - theta[0]);
	
	READ_INT_TYPE resimulation_count = 0;

	//simulating...
	model.startSimulation(&sampler, theta);
	for (READ_INT_TYPE i = 0; i < N; i++) {
		while (!model.simulate(i, read, sid)) { ++resimulation_count; }
		read.write(n_os, os);
		++counts[sid];
		if ((i + 1) % 1000000 == 0 && verbose) cout<<"GEN "<< i + 1<< endl;
	}
	model.finishSimulation();

	cout<< "Total number of resimulation is "<< resimulation_count<< endl;
}

void calcExpressionValues(const vector<double>& theta, const vector<double>& eel, vector<double>& tpm, vector<double>& fpkm) {
	double denom;
	vector<double> frac;

	//calculate fraction of count over all mappabile reads
	denom = 0.0;
	frac.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++)
	  if (eel[i] >= EPSILON) {
	    frac[i] = theta[i];
	    denom += frac[i];
	  }
	general_assert(denom > 0, "No alignable reads?!");
	for (int i = 1; i <= M; i++) frac[i] /= denom;

	//calculate FPKM
	fpkm.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++)
		if (eel[i] >= EPSILON) fpkm[i] = frac[i] * 1e9 / eel[i];

	//calculate TPM
	tpm.assign(M + 1, 0.0);
	denom = 0.0;
	for (int i = 1; i <= M; i++) denom += fpkm[i];
	for (int i = 1; i <= M; i++) tpm[i] = fpkm[i] / denom * 1e6;
}

void writeResFiles(char* outFN) {
	FILE *fo;
	vector<int> tlens;
	vector<double> fpkm, tpm, isopct;
	vector<double> glens, gene_eels, gene_counts, gene_tpm, gene_fpkm;

	for (int i = 1; i <= M; i++)
		general_assert(eel[i] > EPSILON || counts[i] <= EPSILON, "An isoform whose effecitve length < " + ftos(MINEEL, 6) + " got sampled!");

	calcExpressionValues(counts, eel, tpm, fpkm);

	//calculate IsoPct, etc.
	isopct.assign(M + 1, 0.0);
	tlens.assign(M + 1, 0);

	glens.assign(m, 0.0); gene_eels.assign(m, 0.0);
	gene_counts.assign(m, 0.0); gene_tpm.assign(m, 0.0); gene_fpkm.assign(m, 0.0);

	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			const Transcript& transcript = transcripts.getTranscriptAt(j);
			tlens[j] = transcript.getLength();

			glens[i] += tlens[j] * tpm[j];
			gene_eels[i] += eel[j] * tpm[j];
			gene_counts[i] += counts[j];
			gene_tpm[i] += tpm[j];
			gene_fpkm[i] += fpkm[j];
		}

		if (gene_tpm[i] < EPSILON) continue;

		for (int j = b; j < e; j++)
			isopct[j] = tpm[j] / gene_tpm[i];
		glens[i] /= gene_tpm[i];
		gene_eels[i] /= gene_tpm[i];
	}

	//isoform level
	sprintf(isoResF, "%s.sim.isoforms.results", outFN);
	fo = fopen(isoResF, "w");
	fprintf(fo, "transcript_id\tgene_id\tlength\teffective_length\tcount\tTPM\tFPKM\tIsoPct\n");
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		fprintf(fo, "%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str(), tlens[i],
				eel[i], counts[i], tpm[i], fpkm[i], isopct[i] * 1e2);
	}
	fclose(fo);

	//gene level
	sprintf(geneResF, "%s.sim.genes.results", outFN);
	fo = fopen(geneResF, "w");
	fprintf(fo, "gene_id\ttranscript_id(s)\tlength\teffective_length\tcount\tTPM\tFPKM\n");
	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		const string& gene_id = transcripts.getTranscriptAt(b).getGeneID();
		fprintf(fo, "%s\t", gene_id.c_str());
		for (int j = b; j < e; j++) {
			fprintf(fo, "%s%c", transcripts.getTranscriptAt(j).getTranscriptID().c_str(), (j < e - 1 ? ',' : '\t'));
		}
		fprintf(fo, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", glens[i], gene_eels[i], gene_counts[i], gene_tpm[i], gene_fpkm[i]);
	}
	fclose(fo);
}

void releaseOutReadStreams() {
	for (int i = 0; i < n_os; i++) {
		((ofstream*)os[i])->close();
		delete os[i];
	}
}

int main(int argc, char* argv[]) {
	bool quiet = false;
	FILE *fi = NULL;

	if (argc != 7 && argc != 8) {
		printf("Usage: rsem-simulate-reads reference_name estimated_model_file estimated_isoform_results theta0 N output_name [-q]\n\n");
		printf("Parameters:\n\n");
		printf("reference_name: The name of RSEM references, which should be already generated by 'rsem-prepare-reference'\n");
		printf("estimated_model_file: This file describes how the RNA-Seq reads will be sequenced given the expression levels. It determines what kind of reads will be simulated (single-end/paired-end, w/o quality score) and includes parameters for fragment length distribution, read start position distribution, sequencing error models, etc. Normally, this file should be learned from real data using 'rsem-calculate-expression'. The file can be found under the 'sample_name.stat' folder with the name of 'sample_name.model'\n");
		printf("estimated_isoform_results: This file contains expression levels for all isoforms recorded in the reference. It can be learned using 'rsem-calculate-expression' from real data. The corresponding file users want to use is 'sample_name.isoforms.results'. If simulating from user-designed expression profile is desired, start from a learned 'sample_name.isoforms.results' file and only modify the 'TPM' column. The simulator only reads the TPM column. But keeping the file format the same is required.\n");
		printf("theta0: This parameter determines the fraction of reads that are coming from background \"noise\" (instead of from a transcript). It can also be estimated using 'rsem-calculate-expression' from real data. Users can find it as the first value of the third line of the file 'sample_name.stat/sample_name.theta'.\n");
		printf("N: The total number of reads to be simulated. If 'rsem-calculate-expression' is executed on a real data set, the total number of reads can be found as the 4th number of the first line of the file 'sample_name.stat/sample_name.cnt'.\n");
		printf("output_name: Prefix for all output files.\n");
		printf("-q: Set it will stop outputting intermediate information.\n\n");
		printf("Outputs:\n\n");
		printf("output_name.sim.isoforms.results, output_name.sim.genes.results: Expression levels estimated by counting where each simulated read comes from.\n\n");
		printf("output_name.fa if single-end without quality score;\noutput_name.fq if single-end with quality score;\noutput_name_1.fa & output_name_2.fa if paired-end without quality score;\noutput_name_1.fq & output_name_2.fq if paired-end with quality score.\n\n");
		printf("Format of the header line: Each simulated read's header line encodes where it comes from. The header line has the format:\n\n");
		printf("\t{>/@}_rid_dir_sid_pos[_insertL]\n\n");
		printf("{>/@}: Either '>' or '@' must appear. '>' appears if FASTA files are generated and '@' appears if FASTQ files are generated\n");
		printf("rid: Simulated read's index, numbered from 0\n");
		printf("dir: The direction of the simulated read. 0 refers to forward strand ('+') and 1 refers to reverse strand ('-')\n");
		printf("sid: Represent which transcript this read is simulated from. It ranges between 0 and M, where M is the total number of transcripts. If sid=0, the read is simulated from the background noise. Otherwise, the read is simulated from a transcript with index sid. Transcript sid's transcript name can be found in the 'transcript_id' column of the 'sample_name.isoforms.results' file (at line sid + 1, line 1 is for column names)\n");
		printf("pos: The start position of the simulated read in strand dir of transcript sid. It is numbered from 0\n");
		printf("insertL: Only appear for paired-end reads. It gives the insert length of the simulated read.\n\n");
		printf("Example:\n\n");
		printf("Suppose we want to simulate 50 millon single-end reads with quality scores and use the parameters learned from [Example](#example). In addition, we set theta0 as 0.2 and output_name as 'simulated_reads'. The command is:\n\n");
		printf("\trsem-simulate-reads /ref/mouse_125 mmliver_single_quals.stat/mmliver_single_quals.model mmliver_single_quals.isoforms.results 0.2 50000000 simulated_reads\n");
		exit(-1);
	}

	if (argc == 8 && !strcmp(argv[7], "-q")) quiet = true;
	verbose = !quiet;

	//load basic files
	sprintf(refF, "%s.seq", argv[1]);
	refs.loadRefs(refF);
	M = refs.getM();
	sprintf(groupF, "%s.grp", argv[1]);
	gi.load(groupF);
	m = gi.getm();
	sprintf(tiF, "%s.ti", argv[1]);
	transcripts.readFrom(tiF);

	//read model type from modelF
	fi = fopen(argv[2], "r");
	if (fi == NULL) { fprintf(stderr, "Cannot open %s! It may not exist.\n", argv[2]); exit(-1); }
	assert(fscanf(fi, "%d", &model_type) == 1);
	fclose(fi);

	theta.assign(M + 1, 0.0);
	theta[0] = atof(argv[4]);
	N = atoi(argv[5]);

	genOutReadStreams(model_type, argv[6]);

	counts.assign(M + 1, 0.0);

	switch(model_type) {
	case 0: simulate<SingleRead, SingleModel>(argv[2], argv[3]); break;
	case 1: simulate<SingleReadQ, SingleQModel>(argv[2], argv[3]); break;
	case 2: simulate<PairedEndRead, PairedEndModel>(argv[2], argv[3]); break;
	case 3: simulate<PairedEndReadQ, PairedEndQModel>(argv[2], argv[3]); break;
	}

	writeResFiles(argv[6]);
	releaseOutReadStreams();

	return 0;
}
