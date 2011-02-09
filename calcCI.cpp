#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<algorithm>

#include "boost/random.hpp"

#include "utils.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Refs.h"
#include "GroupInfo.h"

using namespace std;

typedef unsigned long bufsize_type;
typedef boost::mt19937 engine_type;
typedef boost::gamma_distribution<> distribution_type;
typedef boost::variate_generator<engine_type&, distribution_type> generator_type;

const int FLOATSIZE = sizeof(float);

struct CIType {
	float lb, ub; // the interval is [lb, ub]

	CIType() { lb = ub = 0.0; }
};

bool quiet;
int model_type;

double confidence;

int nC, cvlen, nSpC, nSamples; // nSpC : number of sample theta vectors per count vector
int fr, to; // each flush, sample fr .. to - 1

int nMB;
bufsize_type size;
float *buffer;
char cvsF[STRLEN], tmpF[STRLEN], command[STRLEN];
ofstream ftmpOut;

int *cvec;
double *theta;
CIType *iso_nrf, *gene_nrf, *iso_tau, *gene_tau;

engine_type engine(time(NULL));
distribution_type **gammas;
generator_type **rgs;

int M, m;
Refs refs;
GroupInfo gi;
char modelF[STRLEN], groupF[STRLEN], refF[STRLEN];

vector<double> eel; //expected effective lengths
double *tau_denoms; // denominators for tau values

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

void flushToTempFile() {
	int gap1 = fr * FLOATSIZE;
	int gap2 = (nSamples - to) * FLOATSIZE;
	float *p = NULL;

	ftmpOut.seekp(0, ios_base::beg);
	for (int i = 0; i < cvlen; i++) {
		p = buffer + i;
		ftmpOut.seekp(gap1, ios_base::cur);
		for (int j = fr; j < to; j++) {
			ftmpOut.write((char*)p, FLOATSIZE);
			p += cvlen;
		}
		ftmpOut.seekp(gap2, ios_base::cur);
	}
}

template<class ModelType>
void sampling() {
	float *p, *ub;
	ifstream fin(cvsF);
	ModelType model;

	model.read(modelF);
	calcExpectedEffectiveLengths<ModelType>(model);

	ftmpOut.open(tmpF, ios::binary);

	fin>>nC>>cvlen;
	assert(cvlen = M + 1);

	nSamples = nC * nSpC;

	fr = to = 0;

	size = bufsize_type(nMB) * 1024 * 1024 / FLOATSIZE / cvlen;
	if (size > (bufsize_type)nSamples) size = nSamples;
	size *= cvlen;
	buffer = new float[size];

	ub = buffer + size;
	p = buffer;

	cvec = new int[cvlen];
	theta = new double[cvlen];
	gammas = new distribution_type*[cvlen];
	rgs = new generator_type*[cvlen];

	tau_denoms = new double[nSamples];
	memset(tau_denoms, 0, sizeof(double) * nSamples);

	double *mw = model.getMW();
	for (int i = 0; i < nC; i++) {
		for (int j = 0; j < cvlen; j++) {
			fin>>cvec[j];
		}

		for (int j = 0; j < cvlen; j++) {
			gammas[j] = new distribution_type(cvec[j]); // need change back before publishing
			rgs[j] = new generator_type(engine, *gammas[j]);
		}

		for (int j = 0; j < nSpC; j++) {
			double sum = 0.0;
			for (int k = 0; k < cvlen; k++) {
				theta[k] = (k == 0 || eel[k] > EPSILON ? (*rgs[k])() : 0.0);
				sum += theta[k];
			}
			assert(sum > 0.0);
			for (int k = 0; k < cvlen; k++) theta[k] /= sum;

			sum = 0.0;
			for (int k = 0; k < cvlen; k++) {
			  theta[k] = (mw[k] < EPSILON ? 0.0 : theta[k] / mw[k]);
			  sum += theta[k];
			}
			assert(sum >= EPSILON);
			for (int k = 0; k < cvlen; k++) theta[k] /= sum;

			*p = (float)theta[0]; ++p;
			assert(1.0 - theta[0] > 0.0);
			for (int k = 1; k < cvlen; k++) {
				if (eel[k] > EPSILON) {
					theta[k] /= (1.0 - theta[0]);
					tau_denoms[to] += theta[k] / eel[k];
				}
				else {
					if (theta[k] != 0.0) { fprintf(stderr, "K=%d Theta_K=%lf\n", k, theta[k]); exit(-1); }
				}

				*p = (float)theta[k];
				++p;
			}
			++to;
			if (p == ub) {
				flushToTempFile();
				p = buffer;
				fr = to;
				if (verbose) { printf("%d vectors are sampled!\n", to); }
			}
		}

		for (int j = 0; j < cvlen; j++) {
			delete gammas[j];
			delete rgs[j];
		}
	}

	if (fr != to) { flushToTempFile(); }

	fin.close();
	ftmpOut.close();

	delete[] buffer;

	delete[] cvec;
	delete[] theta;
	delete[] gammas;
	delete[] rgs;

	if (verbose) { printf("Sampling is finished!\n"); }
}

void calcCI(int nSamples, float *samples, float &lb, float &ub) {
	int p, q; // p pointer for lb, q pointer for ub;
	int newp, newq;
	int threshold = nSamples - (int(confidence * nSamples - 1e-8) + 1);
	int nOutside = 0;

	sort(samples, samples + nSamples);

	p = 0; q = nSamples - 1;
	newq = nSamples - 1;
	do {
		q = newq;
		while (newq > 0 && samples[newq - 1] == samples[newq]) newq--;
		newq--;
	} while (newq >= 0 && nSamples - (newq + 1) <= threshold);

	nOutside = nSamples - (q + 1);

	lb = -1e30; ub = 1e30;
	do {
		if (samples[q] - samples[p] < ub - lb) {
			lb = samples[p];
			ub = samples[q];
		}

		newp = p;
		while (newp < nSamples - 1 && samples[newp] == samples[newp + 1]) newp++;
		newp++;
		if (newp <= threshold) {
			nOutside += newp - p;
			p = newp;
			while (nOutside > threshold && q < nSamples - 1) {
				newq = q + 1;
				while (newq < nSamples - 1 && samples[newq] == samples[newq + 1]) newq++;
				nOutside -= newq - q;
				q = newq;
			}
			assert(nOutside <= threshold);
		}
		else p = newp;
	} while (p <= threshold);
}

void generateResults(char* imdName) {
	float *izsamples, *gzsamples, *itsamples, *gtsamples;
	ifstream fin;
	FILE *fo;
	char outF[STRLEN];

	iso_nrf = new CIType[M + 1];
	iso_tau = new CIType[M + 1];
	gene_nrf = new CIType[m];
	gene_tau = new CIType[m];

	izsamples = new float[nSamples];
	itsamples = new float[nSamples];
	gzsamples = new float[nSamples];
	gtsamples = new float[nSamples];

	fin.open(tmpF, ios::binary);

	for (int k = 0; k < nSamples; k++) fin.read((char*)(&izsamples[k]), FLOATSIZE);
	calcCI(nSamples, izsamples, iso_nrf[0].lb, iso_nrf[0].ub);

	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		memset(gzsamples, 0, FLOATSIZE * nSamples);
		memset(gtsamples, 0, FLOATSIZE * nSamples);
		for (int j = b; j < e; j++) {
			for (int k = 0; k < nSamples; k++) {
				fin.read((char*)(&izsamples[k]), FLOATSIZE);
				if (eel[j] > EPSILON && tau_denoms[k] > EPSILON) { itsamples[k] = izsamples[k] / eel[j] / tau_denoms[k]; }
				else {
					if (izsamples[k] != 0.0) { fprintf(stderr, "K=%d, IZSAMPLES_K=%lf\n", k, izsamples[k]); exit(-1); }
					itsamples[k] = 0.0;
				}
				gzsamples[k] += izsamples[k];
				gtsamples[k] += itsamples[k];
			}
			calcCI(nSamples, izsamples, iso_nrf[j].lb, iso_nrf[j].ub);
			calcCI(nSamples, itsamples, iso_tau[j].lb, iso_tau[j].ub);
		}
		calcCI(nSamples, gzsamples, gene_nrf[i].lb, gene_nrf[i].ub);
		calcCI(nSamples, gtsamples, gene_tau[i].lb, gene_tau[i].ub);

		if (verbose && (i + 1) % 1000 == 0) { printf("%d genes are done!\n", i + 1); }
	}

	fin.close();

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_nrf[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_nrf[i].ub, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_tau[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_tau[i].ub, (i < M ? '\t' : '\n'));
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_nrf[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_nrf[i].ub, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_tau[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_tau[i].ub, (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	printf("CI of noise isoform is [%.6g, %.6g]\n", iso_nrf[0].lb, iso_nrf[0].ub);

	delete[] izsamples;
	delete[] itsamples;
	delete[] gzsamples;
	delete[] gtsamples;

	delete[] iso_nrf;
	delete[] iso_tau;
	delete[] gene_nrf;
	delete[] gene_tau;

	if (verbose) { printf("All credibility intervals are calculated!\n"); }

        sprintf(outF, "%s.tau_denoms", imdName);
        fo = fopen(outF, "w");
        fprintf(fo, "%d\n", nSamples);
        for (int i = 0; i < nSamples; i++) fprintf(fo, "%.15g ", tau_denoms[i]);
        fprintf(fo, "\n");
        fclose(fo);

}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		printf("Usage: rsem-calculate-credibility-intervals reference_name sample_name imdName confidence nSpC nMB[-q]\n");
		exit(-1);
	}

	confidence = atof(argv[4]);
	nSpC = atoi(argv[5]);
	nMB = atoi(argv[6]);

	quiet = false;
	if (argc > 7 && !strcmp(argv[7], "-q")) {
		quiet = true;
	}
	verbose = !quiet;

	sprintf(modelF, "%s.model", argv[2]);
	FILE *fi = fopen(modelF, "r");
	if (fi == NULL) { fprintf(stderr, "Cannot open %s!\n", modelF); exit(-1); }
	fscanf(fi, "%d", &model_type);
	fclose(fi);

	sprintf(refF, "%s.seq", argv[1]);
	refs.loadRefs(refF, 1);
	M = refs.getM();
	sprintf(groupF, "%s.grp", argv[1]);
	gi.load(groupF);
	m = gi.getm();

	sprintf(tmpF, "%s.tmp", argv[3]);
	sprintf(cvsF, "%s.countvectors", argv[3]);

	switch(model_type) {
	case 0 : sampling<SingleModel>(); break;
	case 1 : sampling<SingleQModel>(); break;
	case 2 : sampling<PairedEndModel>(); break;
	case 3 : sampling<PairedEndQModel>(); break;
	}

	generateResults(argv[3]);

	delete[] tau_denoms;

	sprintf(command, "rm -f %s", tmpF);
	int status = system(command);
	if (status != 0) {
	  fprintf(stderr, "Cannot delete %s!\n", tmpF);
	  exit(-1);
	}

	return 0;
}
