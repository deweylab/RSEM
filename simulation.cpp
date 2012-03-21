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

READ_INT_TYPE N;
int model_type, M, m;

Refs refs;
GroupInfo gi;
Transcripts transcripts;

double *theta, *counts;
vector<double> eel;

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

  eel.clear();
  eel.resize(M + 1, 0.0);
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
	double tau;
	double denom = 0.0;
	for (int i = 1; i <= M; i++) {
	  getline(fin, line);
	  size_t pos = 0;
	  for (int j = 0; j < 2; j++) pos = line.find_first_of('\t', pos) + 1;
	  size_t pos2 = line.find_first_of('\t', pos);
	  if (pos2 == string::npos) pos2 = line.length();
	  tau = atof(line.substr(pos, pos2 - pos).c_str());
	  theta[i] = tau * eel[i];
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

void writeResFiles(char* outFN) {
	FILE *fo;
	double denom;

	//calculate tau values
	double *tau = new double[M + 1];
	memset(tau, 0, sizeof(double) * (M + 1));
	denom = 0.0;
	for (int i = 1; i <= M; i++) 
		if (eel[i] > EPSILON) {
			tau[i] = counts[i] / eel[i];
			denom += tau[i];
		}
		else {
		    if (counts[i] > EPSILON) { printf("Warning: An isoform which EEL is less than %.6g gets sampled!\n", MINEEL); }
		}
	assert(denom > 0.0);
	for (int i = 1; i <= M; i++) tau[i] /= denom;

	//isoform level
	sprintf(isoResF, "%s.sim.isoforms.results", outFN);
	fo = fopen(isoResF, "w");
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		fprintf(fo, "%s\t%.2f\t%.15g", transcript.getTranscriptID().c_str(), counts[i], tau[i]);
		
		if (transcript.getLeft() != "") { fprintf(fo, "\t%s", transcript.getLeft().c_str()); }
		fprintf(fo, "\n");
	}
	fclose(fo);

	//gene level
	sprintf(geneResF, "%s.sim.genes.results", outFN);
	fo = fopen(geneResF, "w");
	for (int i = 0; i < m; i++) {
	  double sum_c = 0.0, sum_tau = 0.0;
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			sum_c += counts[j];
			sum_tau += tau[j];
		}
		const string& gene_id = transcripts.getTranscriptAt(b).getGeneID();
		fprintf(fo, "%s\t%.2f\t%.15g\t", gene_id.c_str(), sum_c, sum_tau);
		for (int j = b; j < e; j++) {
			fprintf(fo, "%s%c", transcripts.getTranscriptAt(j).getTranscriptID().c_str(), (j < e - 1 ? ',' : '\n'));
		}
	}
	fclose(fo);

	delete[] tau;
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

	theta = new double[M + 1];
	theta[0] = atof(argv[4]);
	N = atoi(argv[5]);

	genOutReadStreams(model_type, argv[6]);

	counts = new double[M + 1];
	memset(counts, 0, sizeof(double) * (M + 1));

	switch(model_type) {
	case 0: simulate<SingleRead, SingleModel>(argv[2], argv[3]); break;
	case 1: simulate<SingleReadQ, SingleQModel>(argv[2], argv[3]); break;
	case 2: simulate<PairedEndRead, PairedEndModel>(argv[2], argv[3]); break;
	case 3: simulate<PairedEndReadQ, PairedEndQModel>(argv[2], argv[3]); break;
	}

	writeResFiles(argv[6]);
	releaseOutReadStreams();

	delete[] theta;
	delete[] counts;

	return 0;
}
