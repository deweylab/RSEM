#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<algorithm>
#include<vector>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"
#include "sampling.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Refs.h"
#include "GroupInfo.h"
#include "WriteResults.h"

#include "Buffer.h"

using namespace std;

bool verbose = true;

struct Params {
	int no;
	FILE *fi;
	engine_type *engine;
	const double *mw;
};

struct CIParams {
	int no;
	int start_gene_id, end_gene_id;
};

struct CIType {
  float lb, ub; // the interval is [lb, ub]
  float cqv; // coefficient of quartile variation
  
  CIType() { lb = ub = cqv = 0.0; }
};

int model_type;

double pseudoC; // pseudo count, default is 1

int nMB;
double confidence;
int nCV, nSpC, nSamples; // nCV: number of count vectors; nSpC: number of theta vectors sampled per count vector; nSamples: nCV * nSpC
int nThreads;

float *l_bars;

char cvsF[STRLEN], tmpF[STRLEN], command[STRLEN];

CIType *tpm, *fpkm;
CIType *iso_tpm = NULL, *iso_fpkm = NULL;
CIType *gene_tpm, *gene_fpkm;

int M, m;
Refs refs;
GroupInfo gi;
char refName[STRLEN], imdName[STRLEN], statName[STRLEN];
char modelF[STRLEN], groupF[STRLEN], refF[STRLEN];

bool alleleS;
int m_trans;
GroupInfo ta;

vector<double> eel; //expected effective lengths

Buffer *buffer;

bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
int rc;

bool hasSeed;
seedType seed;

CIParams *ciParamsArray;

void* sample_theta_from_c(void* arg) {
	int *cvec;
	double *theta;
	float *tpm;
	gamma_dist **gammas;
	gamma_generator **rgs;

	Params *params = (Params*)arg;
	FILE *fi = params->fi;
	const double *mw = params->mw;

	cvec = new int[M + 1];
	theta = new double[M + 1];
	gammas = new gamma_dist*[M + 1];
	rgs = new gamma_generator*[M + 1];
	tpm = new float[M + 1];
	float l_bar; // the mean transcript length over the sample

	int cnt = 0;
	while (fscanf(fi, "%d", &cvec[0]) == 1) {
		for (int j = 1; j <= M; j++) assert(fscanf(fi, "%d", &cvec[j]) == 1);
		assert(cvec[0] >= 0);

		++cnt;

		for (int j = 0; j <= M; j++) {
		  gammas[j] = NULL; rgs[j] = NULL;
		  if (cvec[j] >= 0) {
		    gammas[j] = new gamma_dist(cvec[j] + pseudoC);
		    rgs[j] = new gamma_generator(*(params->engine), *gammas[j]);
		  }
		}
		  
		for (int i = 0; i < nSpC; i++) {
			double sum = 0.0;
			for (int j = 0; j <= M; j++) {
				theta[j] = ((j == 0 || (cvec[j] >= 0 && eel[j] >= EPSILON && mw[j] >= EPSILON)) ? (*rgs[j])() / mw[j] : 0.0);
				sum += theta[j];
			}
			assert(sum >= EPSILON);
			for (int j = 0; j <= M; j++) theta[j] /= sum;

			sum = 0.0;
			tpm[0] = 0.0;
			for (int j = 1; j <= M; j++)
				if (eel[j] >= EPSILON) {
					tpm[j] = theta[j] / eel[j];
					sum += tpm[j];
				}
				else assert(theta[j] < EPSILON);
			assert(sum >= EPSILON);
			l_bar = 0.0; // store mean effective length of the sample
			for (int j = 1; j <= M; j++) { tpm[j] /= sum; l_bar += tpm[j] * eel[j]; tpm[j] *= 1e6; }
			buffer->write(l_bar, tpm + 1); // ommit the first element in tpm
		}

		for (int j = 0; j <= M; j++) {
		  if (gammas[j] != NULL) delete gammas[j];
		  if (rgs[j] != NULL) delete rgs[j];
		}

		if (verbose && cnt % 100 == 0) { printf("Thread %d, %d count vectors are processed!\n", params->no, cnt); }
	}

	delete[] cvec;
	delete[] theta;
	delete[] gammas;
	delete[] rgs;
	delete[] tpm;

	return NULL;
}

template<class ModelType>
void sample_theta_vectors_from_count_vectors() {
	ModelType model;
	model.read(modelF);
	calcExpectedEffectiveLengths<ModelType>(M, refs, model, eel);

	int num_threads = min(nThreads, nCV);

	buffer = new Buffer(nMB, nSamples, M, l_bars, tmpF);

	paramsArray = new Params[num_threads];
	threads = new pthread_t[num_threads];

	char inpF[STRLEN];
	hasSeed ? engineFactory::init(seed) : engineFactory::init();
	for (int i = 0; i < num_threads; i++) {
		paramsArray[i].no = i;
		sprintf(inpF, "%s%d", cvsF, i);
		paramsArray[i].fi = fopen(inpF, "r");
		paramsArray[i].engine = engineFactory::new_engine();
		paramsArray[i].mw = model.getMW();
	}
	engineFactory::finish();

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (int i = 0; i < num_threads; i++) {
		rc = pthread_create(&threads[i], &attr, &sample_theta_from_c, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) in sample_theta_vectors_from_count_vectors!");
	}
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) in sample_theta_vectors_from_count_vectors!");
	}

	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	for (int i = 0; i < num_threads; i++) {
		fclose(paramsArray[i].fi);
		delete paramsArray[i].engine;
	}
	delete[] paramsArray;

	delete buffer; // Must delete here, force the content left in the buffer be written into the disk

	if (verbose) { printf("Sampling is finished!\n"); }
}

void calcCI(int nSamples, float *samples, CIType& ci) {
	int p, q; // p pointer for lb, q pointer for ub;
	int newp, newq;
	int threshold = nSamples - (int(confidence * nSamples - 1e-8) + 1);
	int nOutside = 0;

	// sort values
	sort(samples, samples + nSamples);

	// calculate credibility interval
	p = 0; q = nSamples - 1;
	newq = nSamples - 1;
	do {
		q = newq;
		while (newq > 0 && samples[newq - 1] == samples[newq]) newq--;
		newq--;
	} while (newq >= 0 && nSamples - (newq + 1) <= threshold);

	nOutside = nSamples - (q + 1);

	ci.lb = -1e30; ci.ub = 1e30;
	do {
		if (samples[q] - samples[p] < ci.ub - ci.lb) {
			ci.lb = samples[p];
			ci.ub = samples[q];
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

	
	// calculate coefficient of quartile variation
	float Q1, Q3; // the first and third quartiles

	// calculate Tukey's hinges
	int quotient = nSamples / 4;
	int residue = nSamples % 4;

	if (residue == 0) {
	  Q1 = (samples[quotient - 1] + samples[quotient]) / 2.0;
	  Q3 = (samples[3 * quotient - 1] + samples[3 * quotient]) / 2.0;
	}
	else if (residue == 3) {
	  Q1 = (samples[quotient] + samples[quotient + 1]) / 2.0;
	  Q3 = (samples[quotient * 3 + 1] + samples[quotient * 3 + 2]) / 2.0;
	}
	else {
	  Q1 = samples[quotient];
	  Q3 = samples[3 * quotient];
	}

	ci.cqv = (Q3 - Q1 > 0.0 ? (Q3 - Q1) / (Q3 + Q1) : 0.0);
}

void* calcCI_batch(void* arg) {
	float *tsamples, *fsamples;
	float *itsamples = NULL, *ifsamples = NULL, *gtsamples, *gfsamples;
	ifstream fin;
	CIParams *ciParams = (CIParams*)arg;
	int curtid, curaid, tid;

	tsamples = new float[nSamples];
	fsamples = new float[nSamples];
	if (alleleS) { 
	  itsamples = new float[nSamples];
	  ifsamples = new float[nSamples];
	}
	gtsamples = new float[nSamples];
	gfsamples = new float[nSamples];

	fin.open(tmpF, ios::binary);
	// minus 1 here for that theta0 is not written!
	streampos pos = streampos(gi.spAt(ciParams->start_gene_id) - 1) * nSamples * FLOATSIZE;
	fin.seekg(pos, ios::beg);

	int cnt = 0;
	if (alleleS) {
	  curtid = curaid = -1;
	  memset(itsamples, 0, FLOATSIZE * nSamples);
	  memset(ifsamples, 0, FLOATSIZE * nSamples);
	}
	for (int i = ciParams->start_gene_id; i < ciParams->end_gene_id; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		memset(gtsamples, 0, FLOATSIZE * nSamples);
		memset(gfsamples, 0, FLOATSIZE * nSamples);
		for (int j = b; j < e; j++) {
			if (alleleS) {
			  tid = ta.gidAt(j);
			  if (curtid != tid) {
			    if (curtid >= 0) {
			      if (j - curaid > 1) {
				calcCI(nSamples, itsamples, iso_tpm[curtid]);
				calcCI(nSamples, ifsamples, iso_fpkm[curtid]);
			      }
			      else {
				iso_tpm[curtid] = tpm[curaid];
				iso_fpkm[curtid] = fpkm[curaid];
			      }
			    }
			    curtid = tid;
			    curaid = j;
			  }
			}

			for (int k = 0; k < nSamples; k++) {
				fin.read((char*)(&tsamples[k]), FLOATSIZE);
				fsamples[k] = 1e3 / l_bars[k] * tsamples[k];
				if (alleleS) {
				  itsamples[k] += tsamples[k];
				  ifsamples[k] += fsamples[k];
				}
				gtsamples[k] += tsamples[k];
				gfsamples[k] += fsamples[k];
			}
			calcCI(nSamples, tsamples, tpm[j]);
			calcCI(nSamples, fsamples, fpkm[j]);
		}

		if (e - b > 1) {
		  calcCI(nSamples, gtsamples, gene_tpm[i]);
		  calcCI(nSamples, gfsamples, gene_fpkm[i]);
		}
		else {
			gene_tpm[i] = tpm[b];
			gene_fpkm[i] = fpkm[b];
		}

		++cnt;
		if (verbose && cnt % 1000 == 0) { printf("In thread %d, %d genes are processed for CI calculation!\n", ciParams->no, cnt); }
	}
	fin.close();

	if (alleleS && (curtid >= 0)) {
	  if (gi.spAt(ciParams->end_gene_id) - curaid > 1) {
	    calcCI(nSamples, itsamples, iso_tpm[curtid]);
	    calcCI(nSamples, ifsamples, iso_fpkm[curtid]);
	  }
	  else {
	    iso_tpm[curtid] = tpm[curaid];
	    iso_fpkm[curtid] = fpkm[curaid];
	  }
	}
	
	delete[] tsamples;
	delete[] fsamples;
	if (alleleS) {
	  delete[] itsamples;
	  delete[] ifsamples;
	}
	delete[] gtsamples;
	delete[] gfsamples;

	return NULL;
}

void calculate_credibility_intervals(char* imdName) {
	FILE *fo;
	char outF[STRLEN];
	int num_threads = nThreads;

	tpm = new CIType[M + 1];
	fpkm = new CIType[M + 1];
	if (alleleS) {
	  iso_tpm = new CIType[m_trans];
	  iso_fpkm = new CIType[m_trans];
	}
	gene_tpm = new CIType[m];
	gene_fpkm = new CIType[m];

	assert(M > 0);
	int quotient = M / num_threads;
	if (quotient < 1) { num_threads = M; quotient = 1; }
	int cur_gene_id = 0;
	int num_isoforms = 0;

	// A just so so strategy for paralleling
	ciParamsArray = new CIParams[num_threads];
	for (int i = 0; i < num_threads; i++) {
		ciParamsArray[i].no = i;
		ciParamsArray[i].start_gene_id = cur_gene_id;
		num_isoforms = 0;

		while ((m - cur_gene_id > num_threads - i - 1) && (i == num_threads - 1 || num_isoforms < quotient)) {
			num_isoforms += gi.spAt(cur_gene_id + 1) - gi.spAt(cur_gene_id);
			++cur_gene_id;
		}

		ciParamsArray[i].end_gene_id = cur_gene_id;
	}

	threads = new pthread_t[num_threads];

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// paralleling
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_create(&threads[i], &attr, &calcCI_batch, (void*)(&ciParamsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) in calculate_credibility_intervals!");
	}
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) in calculate_credibility_intervals!");
	}

	// releasing resources

	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	delete[] ciParamsArray;

	alleleS ? sprintf(outF, "%s.allele_res", imdName) : sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", tpm[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", tpm[i].ub, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", tpm[i].cqv, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", fpkm[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", fpkm[i].ub, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
	  fprintf(fo, "%.6g%c", fpkm[i].cqv, (i < M ? '\t' : '\n'));
	fclose(fo);

	if (alleleS) {
	  //isoform level results
	  sprintf(outF, "%s.iso_res", imdName);
	  fo = fopen(outF, "a");
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_tpm[i].lb, (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_tpm[i].ub, (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_tpm[i].cqv, (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_fpkm[i].lb, (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_fpkm[i].ub, (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.6g%c", iso_fpkm[i].cqv, (i < m_trans - 1 ? '\t' : '\n'));
	  fclose(fo);
	}

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_tpm[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_tpm[i].ub, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_tpm[i].cqv, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_fpkm[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_fpkm[i].ub, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.6g%c", gene_fpkm[i].cqv, (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	delete[] tpm;
	delete[] fpkm;
	if (alleleS) { 
	  delete[] iso_tpm; 
	  delete[] iso_fpkm; 
	}
	delete[] gene_tpm;
	delete[] gene_fpkm;

	if (verbose) { printf("All credibility intervals are calculated!\n"); }
}

int main(int argc, char* argv[]) {
	if (argc < 8) {
		printf("Usage: rsem-calculate-credibility-intervals reference_name imdName statName confidence nCV nSpC nMB [-p #Threads] [--seed seed] [--pseudo-count pseudo_count] [-q]\n");
		exit(-1);
	}

	strcpy(refName, argv[1]);
	strcpy(imdName, argv[2]);
	strcpy(statName, argv[3]);

	confidence = atof(argv[4]);
	nCV = atoi(argv[5]);
	nSpC = atoi(argv[6]);
	nMB = atoi(argv[7]);

	nThreads = 1;
	quiet = false;
	hasSeed = false;
	pseudoC = 1.0;
	for (int i = 8; i < argc; i++) {
		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--seed")) {
		  hasSeed = true;
		  int len = strlen(argv[i + 1]);
		  seed = 0;
		  for (int k = 0; k < len; k++) seed = seed * 10 + (argv[i + 1][k] - '0');
		}
		if (!strcmp(argv[i], "--pseudo-count")) pseudoC = atof(argv[i + 1]);
		if (!strcmp(argv[i], "-q")) quiet = true;
	}
	verbose = !quiet;

	sprintf(refF, "%s.seq", refName);
	refs.loadRefs(refF, 1);
	M = refs.getM();

	sprintf(groupF, "%s.grp", refName);
	gi.load(groupF);
	m = gi.getm();

	// allele-specific 
	alleleS = isAlleleSpecific(refName, NULL, &ta);
	if (alleleS) m_trans = ta.getm();

	nSamples = nCV * nSpC;
	assert(nSamples > 0 && M > 0); // for Buffter.h: (bufsize_type)nSamples
	l_bars = new float[nSamples];

	sprintf(tmpF, "%s.tmp", imdName);
	sprintf(cvsF, "%s.countvectors", imdName);

	sprintf(modelF, "%s.model", statName);
	FILE *fi = fopen(modelF, "r");
	general_assert(fi != NULL, "Cannot open " + cstrtos(modelF) + "!");
	assert(fscanf(fi, "%d", &model_type) == 1);
	fclose(fi);

	// Phase I
	switch(model_type) {
	case 0 : sample_theta_vectors_from_count_vectors<SingleModel>(); break;
	case 1 : sample_theta_vectors_from_count_vectors<SingleQModel>(); break;
	case 2 : sample_theta_vectors_from_count_vectors<PairedEndModel>(); break;
	case 3 : sample_theta_vectors_from_count_vectors<PairedEndQModel>(); break;
	}

	// Phase II
	calculate_credibility_intervals(imdName);

	delete l_bars;

	return 0;
}
