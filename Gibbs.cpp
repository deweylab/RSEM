#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<sstream>
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

using namespace std;

struct Params {
	int no, nsamples;
	FILE *fo;
	engine_type *engine;
	double *pme_c, *pve_c; //posterior mean and variance vectors on counts
	double *pme_theta;
};


struct Item {
	int sid;
	double conprb;

	Item(int sid, double conprb) {
		this->sid = sid;
		this->conprb = conprb;
	}
};

int nThreads;

int model_type;
int m, M;
READ_INT_TYPE N0, N1;
HIT_INT_TYPE nHits;
double totc;
int BURNIN, NSAMPLES, GAP;
char imdName[STRLEN], statName[STRLEN];
char thetaF[STRLEN], ofgF[STRLEN], groupF[STRLEN], refF[STRLEN], modelF[STRLEN];
char cvsF[STRLEN];

Refs refs;
GroupInfo gi;

vector<HIT_INT_TYPE> s;
vector<Item> hits;

vector<double> theta;

vector<double> pme_c, pve_c; //global posterior mean and variance vectors on counts
vector<double> pme_theta, eel;

bool var_opt;
bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
void *status;
int rc;

void load_data(char* reference_name, char* statName, char* imdName) {
	ifstream fin;
	string line;
	int tmpVal;

	//load reference file
	sprintf(refF, "%s.seq", reference_name);
	refs.loadRefs(refF, 1);
	M = refs.getM();

	//load groupF
	sprintf(groupF, "%s.grp", reference_name);
	gi.load(groupF);
	m = gi.getm();

	//load thetaF
	sprintf(thetaF, "%s.theta",statName);
	fin.open(thetaF);
	general_assert(fin.is_open(), "Cannot open " + cstrtos(thetaF) + "!");
	fin>>tmpVal;
	general_assert(tmpVal == M + 1, "Number of transcripts is not consistent in " + cstrtos(refF) + " and " + cstrtos(thetaF) + "!");
	theta.assign(M + 1, 0);
	for (int i = 0; i <= M; i++) fin>>theta[i];
	fin.close();

	//load ofgF;
	sprintf(ofgF, "%s.ofg", imdName);
	fin.open(ofgF);
	general_assert(fin.is_open(), "Cannot open " + cstrtos(ofgF) + "!");
	fin>>tmpVal>>N0;
	general_assert(tmpVal == M, "M in " + cstrtos(ofgF) + " is not consistent with " + cstrtos(refF) + "!");
	getline(fin, line);

	s.clear(); hits.clear();
	s.push_back(0);
	while (getline(fin, line)) {
		istringstream strin(line);
		int sid;
		double conprb;

		while (strin>>sid>>conprb) {
			hits.push_back(Item(sid, conprb));
		}
		s.push_back(hits.size());
	}
	fin.close();

	N1 = s.size() - 1;
	nHits = hits.size();

	totc = N0 + N1 + (M + 1);

	if (verbose) { printf("Loading Data is finished!\n"); }
}

// assign threads
void init() {
	int quotient, left;
	char outF[STRLEN];

	quotient = NSAMPLES / nThreads;
	left = NSAMPLES % nThreads;

	sprintf(cvsF, "%s.countvectors", imdName);
	paramsArray = new Params[nThreads];
	threads = new pthread_t[nThreads];

	for (int i = 0; i < nThreads; i++) {
		paramsArray[i].no = i;

		paramsArray[i].nsamples = quotient;
		if (i < left) paramsArray[i].nsamples++;

		sprintf(outF, "%s%d", cvsF, i);
		paramsArray[i].fo = fopen(outF, "w");

		paramsArray[i].engine = engineFactory::new_engine();
		paramsArray[i].pme_c = new double[M + 1];
		memset(paramsArray[i].pme_c, 0, sizeof(double) * (M + 1));
		paramsArray[i].pve_c = new double[M + 1];
		memset(paramsArray[i].pve_c, 0, sizeof(double) * (M + 1));
		paramsArray[i].pme_theta = new double[M + 1];
		memset(paramsArray[i].pme_theta, 0, sizeof(double) * (M + 1));
	}

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	if (verbose) { printf("Initialization finished!\n"); }
}

//sample theta from Dir(1)
void sampleTheta(engine_type& engine, vector<double>& theta) {
	gamma_dist gm(1);
	gamma_generator gmg(engine, gm);
	double denom;

	theta.assign(M + 1, 0);
	denom = 0.0;
	for (int i = 0; i <= M; i++) {
		theta[i] = gmg();
		denom += theta[i];
	}
	assert(denom > EPSILON);
	for (int i = 0; i <= M; i++) theta[i] /= denom;
}

void writeCountVector(FILE* fo, vector<int>& counts) {
	for (int i = 0; i < M; i++) {
		fprintf(fo, "%d ", counts[i]);
	}
	fprintf(fo, "%d\n", counts[M]);
}

void* Gibbs(void* arg) {
	int CHAINLEN;
	HIT_INT_TYPE len, fr, to;
	Params *params = (Params*)arg;

	vector<double> theta;
	vector<int> z, counts;
	vector<double> arr;

	uniform01 rg(*params->engine);

	// generate initial state
	sampleTheta(*params->engine, theta);

	z.assign(N1, 0);

	counts.assign(M + 1, 1); // 1 pseudo count
	counts[0] += N0;

	for (READ_INT_TYPE i = 0; i < N1; i++) {
		fr = s[i]; to = s[i + 1];
		len = to - fr;
		arr.assign(len, 0);
		for (HIT_INT_TYPE j = fr; j < to; j++) {
			arr[j - fr] = theta[hits[j].sid] * hits[j].conprb;
			if (j > fr) arr[j - fr] += arr[j - fr - 1];  // cumulative
		}
		z[i] = hits[fr + sample(rg, arr, len)].sid;
		++counts[z[i]];
	}

	// Gibbs sampling
	CHAINLEN = 1 + (params->nsamples - 1) * GAP;
	for (int ROUND = 1; ROUND <= BURNIN + CHAINLEN; ROUND++) {

		for (READ_INT_TYPE i = 0; i < N1; i++) {
			--counts[z[i]];
			fr = s[i]; to = s[i + 1]; len = to - fr;
			arr.assign(len, 0);
			for (HIT_INT_TYPE j = fr; j < to; j++) {
				arr[j - fr] = counts[hits[j].sid] * hits[j].conprb;
				if (j > fr) arr[j - fr] += arr[j - fr - 1]; //cumulative
			}
			z[i] = hits[fr + sample(rg, arr, len)].sid;
			++counts[z[i]];
		}

		if (ROUND > BURNIN) {
			if ((ROUND - BURNIN - 1) % GAP == 0) {
				writeCountVector(params->fo, counts);
				for (int i = 0; i <= M; i++) {
					params->pme_c[i] += counts[i] - 1;
					params->pve_c[i] += (counts[i] - 1) * (counts[i] - 1);
					params->pme_theta[i] += counts[i] / totc;
				}
			}
		}

		if (verbose && ROUND % 100 == 0) { printf("Thread %d, ROUND %d is finished!\n", params->no, ROUND); }
	}

	return NULL;
}

void release() {
//	char inpF[STRLEN], command[STRLEN];
	string line;

	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	pme_c.assign(M + 1, 0);
	pve_c.assign(M + 1, 0);
	pme_theta.assign(M + 1, 0);
	for (int i = 0; i < nThreads; i++) {
		fclose(paramsArray[i].fo);
		delete paramsArray[i].engine;
		for (int j = 0; j <= M; j++) {
			pme_c[j] += paramsArray[i].pme_c[j];
			pve_c[j] += paramsArray[i].pve_c[j];
			pme_theta[j] += paramsArray[i].pme_theta[j];
		}
		delete[] paramsArray[i].pme_c;
		delete[] paramsArray[i].pve_c;
		delete[] paramsArray[i].pme_theta;
	}
	delete[] paramsArray;


	for (int i = 0; i <= M; i++) {
		pme_c[i] /= NSAMPLES;
		pve_c[i] = (pve_c[i] - NSAMPLES * pme_c[i] * pme_c[i]) / (NSAMPLES - 1);
		pme_theta[i] /= NSAMPLES;
	}
}

template<class ModelType>
void calcExpectedEffectiveLengths(ModelType& model) {
	int lb, ub, span;
	double *pdf = NULL, *cdf = NULL, *clen = NULL; // clen[i] = \sigma_{j=1}^{i}pdf[i]*(lb+i)
  
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

template<class ModelType>
void writeEstimatedParameters(char* modelF, char* imdName) {
	ModelType model;
	double denom;
	char outF[STRLEN];
	FILE *fo;

	model.read(modelF);

	calcExpectedEffectiveLengths<ModelType>(model);

	denom = pme_theta[0];
	for (int i = 1; i <= M; i++)
	  if (eel[i] < EPSILON) pme_theta[i] = 0.0;
	  else denom += pme_theta[i];

	general_assert(denom >= EPSILON, "No Expected Effective Length is no less than " + ftos(MINEEL, 6) + "?!");

	for (int i = 0; i <= M; i++) pme_theta[i] /= denom;

	denom = 0.0;
	double *mw = model.getMW();
	for (int i = 0; i <= M; i++) {
	  pme_theta[i] = (mw[i] < EPSILON ? 0.0 : pme_theta[i] / mw[i]);
	  denom += pme_theta[i];
	}
	assert(denom >= EPSILON);
	for (int i = 0; i <= M; i++) pme_theta[i] /= denom;

	//calculate tau values
	double *tau = new double[M + 1];
	memset(tau, 0, sizeof(double) * (M + 1));

	denom = 0.0;
	for (int i = 1; i <= M; i++) 
	  if (eel[i] > EPSILON) {
	    tau[i] = pme_theta[i] / eel[i];
	    denom += tau[i];
	  }

	general_assert(denom >= EPSILON, "No alignable reads?!");

	for (int i = 1; i <= M; i++) {
		tau[i] /= denom;
	}

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");

	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", pme_c[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.15g%c", tau[i], (i < M ? '\t' : '\n'));

	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");

	for (int i = 0; i < m; i++) {
		double sumC = 0.0; //  sum of pme counts
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			sumC += pme_c[j];
		}
		fprintf(fo, "%.15g%c", sumC, (i < m - 1 ? '\t' : '\n'));
	}
	for (int i = 0; i < m; i++) {
		double sumT = 0.0; //  sum of tau values
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			sumT += tau[j];
		}
		fprintf(fo, "%.15g%c", sumT, (i < m - 1 ? '\t' : '\n'));
	}
	fclose(fo);

	delete[] tau;

	if (verbose) { printf("Gibbs based expression values are written!\n"); }
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		printf("Usage: rsem-run-gibbs-multi reference_name sample_name sampleToken BURNIN NSAMPLES GAP [-p #Threads] [--var] [-q]\n");
		exit(-1);
	}

	BURNIN = atoi(argv[4]);
	NSAMPLES = atoi(argv[5]);
	GAP = atoi(argv[6]);
	sprintf(imdName, "%s.temp/%s", argv[2], argv[3]);
	sprintf(statName, "%s.stat/%s", argv[2], argv[3]);
	load_data(argv[1], statName, imdName);

	nThreads = 1;
	var_opt = false;
	quiet = false;

	for (int i = 7; i < argc; i++) {
		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--var")) var_opt = true;
		if (!strcmp(argv[i], "-q")) quiet = true;
	}
	verbose = !quiet;

	assert(NSAMPLES > 1); // Otherwise, we cannot calculate posterior variance

	if (nThreads > NSAMPLES) {
		nThreads = NSAMPLES;
		printf("Warning: Number of samples is less than number of threads! Change the number of threads to %d!\n", nThreads);
	}

	if (verbose) printf("Gibbs started!\n");

	init();
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_create(&threads[i], &attr, Gibbs, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0)!");
	}
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_join(threads[i], &status);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0)!");
	}
	release();

	if (verbose) printf("Gibbs finished!\n");

	sprintf(modelF, "%s.model", statName);
	FILE *fi = fopen(modelF, "r");
	general_assert(fi != NULL, "Cannot open " + cstrtos(modelF) + "!");
	assert(fscanf(fi, "%d", &model_type) == 1);
	fclose(fi);

	switch(model_type) {
	case 0 : writeEstimatedParameters<SingleModel>(modelF, imdName); break;
	case 1 : writeEstimatedParameters<SingleQModel>(modelF, imdName); break;
	case 2 : writeEstimatedParameters<PairedEndModel>(modelF, imdName); break;
	case 3 : writeEstimatedParameters<PairedEndQModel>(modelF, imdName); break;
	}

	if (var_opt) {
		char varF[STRLEN];

		sprintf(varF, "%s.var", statName);
		FILE *fo = fopen(varF, "w");
		general_assert(fo != NULL, "Cannot open " + cstrtos(varF) + "!");
		for (int i = 0; i < m; i++) {
			int b = gi.spAt(i), e = gi.spAt(i + 1), number_of_isoforms = e - b;
			for (int j = b; j < e; j++) {
				fprintf(fo, "%s\t%d\t%.15g\t%.15g\n", refs.getRef(j).getName().c_str(), number_of_isoforms, pme_c[j], pve_c[j]);
			}
		}
		fclose(fo);
	}

	return 0;
}
