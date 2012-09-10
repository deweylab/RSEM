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
	  theta[i] = tpm * eel[i];
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
		printf("Usage: rsem-simulate-reads reference_name estimated_model_file estimated_isoform_results theta0 N output_name [-q]\n");
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
