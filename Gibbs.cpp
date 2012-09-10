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
  double *pme_tpm, *pme_fpkm;
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

vector<double> eel;
double *mw;

vector<double> pme_c, pve_c; //global posterior mean and variance vectors on counts
vector<double> pme_tpm, pme_fpkm;

bool var_opt;
bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
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

template<class ModelType>
void init_model_related(char* modelF) {
	ModelType model;
	model.read(modelF);

	calcExpectedEffectiveLengths<ModelType>(model);
	memcpy(mw, model.getMW(), sizeof(double) * (M + 1)); // otherwise, after exiting this procedure, mw becomes undefined
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
		paramsArray[i].pme_tpm = new double[M + 1];
		memset(paramsArray[i].pme_tpm, 0, sizeof(double) * (M + 1));
		paramsArray[i].pme_fpkm = new double[M + 1];
		memset(paramsArray[i].pme_fpkm, 0, sizeof(double) * (M + 1));
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

void polishTheta(vector<double>& theta, const vector<double>& eel, const double* mw) {
	double sum = 0.0;

	/* The reason that for noise gene, mw value is 1 is :
	 * currently, all masked positions are for poly(A) sites, which in theory should be filtered out.
	 * So the theta0 does not containing reads from any masked position
	 */

	for (int i = 0; i <= M; i++) {
		// i == 0, mw[i] == 1
		if (i > 0 && (mw[i] < EPSILON || eel[i] < EPSILON)) {
			theta[i] = 0.0;
			continue;
		}
		theta[i] = theta[i] / mw[i];
		sum += theta[i];
	}
	// currently is OK, since no transcript should be masked totally, only the poly(A) tail related part will be masked
	general_assert(sum >= EPSILON, "No effective length is no less than" + ftos(MINEEL, 6) + " !");
	for (int i = 0; i <= M; i++) theta[i] /= sum;
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

void* Gibbs(void* arg) {
	int CHAINLEN;
	HIT_INT_TYPE len, fr, to;
	Params *params = (Params*)arg;

	vector<double> theta, tpm, fpkm;
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
				for (int i = 0; i <= M; i++) theta[i] = counts[i] / totc;
				polishTheta(theta, eel, mw);
				calcExpressionValues(theta, eel, tpm, fpkm);
				for (int i = 0; i <= M; i++) {
					params->pme_c[i] += counts[i] - 1;
					params->pve_c[i] += (counts[i] - 1) * (counts[i] - 1);
					params->pme_tpm[i] += tpm[i];
					params->pme_fpkm[i] += fpkm[i];
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
	pme_tpm.assign(M + 1, 0);
	pme_fpkm.assign(M + 1, 0);
	for (int i = 0; i < nThreads; i++) {
		fclose(paramsArray[i].fo);
		delete paramsArray[i].engine;
		for (int j = 0; j <= M; j++) {
			pme_c[j] += paramsArray[i].pme_c[j];
			pve_c[j] += paramsArray[i].pve_c[j];
			pme_tpm[j] += paramsArray[i].pme_tpm[j];
			pme_fpkm[j] += paramsArray[i].pme_fpkm[j];
		}
		delete[] paramsArray[i].pme_c;
		delete[] paramsArray[i].pve_c;
		delete[] paramsArray[i].pme_tpm;
		delete[] paramsArray[i].pme_fpkm;
	}
	delete[] paramsArray;


	for (int i = 0; i <= M; i++) {
		pme_c[i] /= NSAMPLES;
		pve_c[i] = (pve_c[i] - NSAMPLES * pme_c[i] * pme_c[i]) / (NSAMPLES - 1);
		pme_tpm[i] /= NSAMPLES;
		pme_fpkm[i] /= NSAMPLES;
	}
}

void writeResults(char* imdName) {
	char outF[STRLEN];
	FILE *fo;

	vector<double> isopct;
	vector<double> gene_counts, gene_tpm, gene_fpkm;

	//calculate IsoPct, etc.
	isopct.assign(M + 1, 0.0);
	gene_counts.assign(m, 0.0); gene_tpm.assign(m, 0.0); gene_fpkm.assign(m, 0.0);

	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			gene_counts[i] += pme_c[j];
			gene_tpm[i] += pme_tpm[j];
			gene_fpkm[i] += pme_fpkm[j];
		}
		if (gene_tpm[i] < EPSILON) continue;
		for (int j = b; j < e; j++)
			isopct[j] = pme_tpm[j] / gene_tpm[i];
	}

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");

	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", pme_c[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", pme_tpm[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", pme_fpkm[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");

	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_counts[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_tpm[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_fpkm[i], (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	if (verbose) { printf("Gibbs based expression values are written!\n"); }
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		printf("Usage: rsem-run-gibbs reference_name imdName statName BURNIN NSAMPLES GAP [-p #Threads] [--var] [-q]\n");
		exit(-1);
	}

	strcpy(imdName, argv[2]);
	strcpy(statName, argv[3]);

	BURNIN = atoi(argv[4]);
	NSAMPLES = atoi(argv[5]);
	GAP = atoi(argv[6]);

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

	load_data(argv[1], statName, imdName);

	sprintf(modelF, "%s.model", statName);
	FILE *fi = fopen(modelF, "r");
	general_assert(fi != NULL, "Cannot open " + cstrtos(modelF) + "!");
	assert(fscanf(fi, "%d", &model_type) == 1);
	fclose(fi);

	mw = new double[M + 1]; // make an extra copy

	switch(model_type) {
	case 0 : init_model_related<SingleModel>(modelF); break;
	case 1 : init_model_related<SingleQModel>(modelF); break;
	case 2 : init_model_related<PairedEndModel>(modelF); break;
	case 3 : init_model_related<PairedEndQModel>(modelF); break;
	}

	if (verbose) printf("Gibbs started!\n");

	init();
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_create(&threads[i], &attr, Gibbs, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0)!");
	}
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0)!");
	}
	release();

	if (verbose) printf("Gibbs finished!\n");
	
	writeResults(imdName);

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

	delete mw; // delete the copy

	return 0;
}
