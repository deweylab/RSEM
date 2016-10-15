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
#include "WriteResults.h"

using namespace std;

bool verbose = true;

struct Params {
  int no, nsamples;
  FILE *fo;
  engine_type *engine;
  double *pme_c, *pve_c; //posterior mean and variance vectors on counts
  double *pme_tpm, *pme_fpkm;
  
  double *pve_c_genes, *pve_c_trans;
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
int M;
READ_INT_TYPE N0, N1;
HIT_INT_TYPE nHits;
double totc;
int BURNIN, NSAMPLES, GAP;
char refName[STRLEN], imdName[STRLEN], statName[STRLEN];
char thetaF[STRLEN], ofgF[STRLEN], refF[STRLEN], modelF[STRLEN];
char cvsF[STRLEN];

Refs refs;

vector<HIT_INT_TYPE> s;
vector<Item> hits;

vector<double> eel;
double *mw;

vector<int> init_counts;
double pseudoC;

vector<double> pme_c, pve_c; //global posterior mean and variance vectors on counts
vector<double> pme_tpm, pme_fpkm;

bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
int rc;

bool hasSeed;
seedType seed;

int m;
char groupF[STRLEN];
GroupInfo gi;

bool alleleS;
int m_trans;
GroupInfo gt, ta;
vector<double> pve_c_genes, pve_c_trans;

// pliu
// if has prior file; file's name; and a vector to save prior parameters
bool has_prior;
char fprior[STRLEN];
vector<double> pseudo_counts;
//////

void load_data(char* refName, char* statName, char* imdName) {
	ifstream fin;
	string line;
	int tmpVal;

	//load reference file
	sprintf(refF, "%s.seq", refName);
	refs.loadRefs(refF, 1);
	M = refs.getM();

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

	if (verbose) { printf("Loading data is finished!\n"); }
}

void load_group_info(char* refName) {
  // Load group info
  sprintf(groupF, "%s.grp", refName);
  gi.load(groupF);
  m = gi.getm();
  
  alleleS = isAlleleSpecific(refName, &gt, &ta); // if allele-specific 
  m_trans = (alleleS ? ta.getm() : 0);

  if (verbose) { printf("Loading group information is finished!\n"); }
}

// Load imdName.omit and initialize the init count vector.
void load_omit_info(const char* imdName) {
  char omitF[STRLEN];
  FILE *fi = NULL;
  int tid;
  
  sprintf(omitF, "%s.omit", imdName);
  fi = fopen(omitF, "r");
  init_counts.assign(M + 1, 0);
  totc = M + 1;
  while (fscanf(fi, "%d", &tid) == 1) {
    init_counts[tid] = -1;
    --totc;
  }
  fclose(fi);
  totc = totc * pseudoC + N0 + N1;
}

// pliu
// load isoform's prior information and re-calculate totc
void load_prior_info(const char* fprior){
	pseudo_counts.assign(M+1, 0.0);
  ifstream fin;
  string line;
  fin.open(fprior);
  for(int i=1; i<=M; ++i){
    double prior;
    getline(fin, line);
    sscanf(line.c_str(), "%lf%*s", &prior);
    if ( init_counts[i] == 0 ){ // not to-be-omitted
      pseudo_counts[i] = prior;
    } 
  }
  fin.close();

	// re-calculate 'totc' by considering prior parameters
	totc = 1;
	for ( int i=1; i<=M; ++i ) {
		if ( init_counts[i] == 0 ) { // not to-be-omitted
			totc += pseudo_counts[i];
		}
	}
	totc += N0 + N1;
}
//////

template<class ModelType>
void init_model_related(char* modelF) {
	ModelType model;
	model.read(modelF);

	calcExpectedEffectiveLengths<ModelType>(M, refs, model, eel);
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

	hasSeed ? engineFactory::init(seed) : engineFactory::init();
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

		paramsArray[i].pve_c_genes = new double[m];
		memset(paramsArray[i].pve_c_genes, 0, sizeof(double) * m);
		
		paramsArray[i].pve_c_trans = NULL;
		if (alleleS) {
		  paramsArray[i].pve_c_trans = new double[m_trans];
		  memset(paramsArray[i].pve_c_trans, 0, sizeof(double) * m_trans);
		}
	}
	engineFactory::finish();

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	if (verbose) { printf("Initialization finished!\n"); }
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

	vector<double> theta, tpm, fpkm;
	vector<int> z, counts(init_counts);
	vector<double> arr;

	uniform_01_generator rg(*params->engine, uniform_01_dist());

	// generate initial state
	theta.assign(M + 1, 0.0);
	z.assign(N1, 0);
	counts[0] += N0;

	for (READ_INT_TYPE i = 0; i < N1; i++) {
		fr = s[i]; to = s[i + 1];
		len = to - fr;
		arr.assign(len, 0);
		for (HIT_INT_TYPE j = fr; j < to; j++) {
			arr[j - fr] = hits[j].conprb;
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
				if ( has_prior ) {
			  	arr[j - fr] = (counts[hits[j].sid] + pseudo_counts[hits[j].sid]) * hits[j].conprb;
				} else {
			  	arr[j - fr] = (counts[hits[j].sid] + pseudoC) * hits[j].conprb;
				}
			  if (j > fr) arr[j - fr] += arr[j - fr - 1]; //cumulative
			}
			z[i] = hits[fr + sample(rg, arr, len)].sid;
			++counts[z[i]];
		}

		if (ROUND > BURNIN) {
			if ((ROUND - BURNIN - 1) % GAP == 0) {
				writeCountVector(params->fo, counts);
				for (int i = 0; i <= M; i++) {
					if ( has_prior ) {
						theta[i] = (counts[i] < 0 ? 0.0 : (counts[i] + pseudo_counts[i]) / totc);
					} else {
						theta[i] = (counts[i] < 0 ? 0.0 : (counts[i] + pseudoC) / totc);
					}
				}
				polishTheta(M, theta, eel, mw);
				calcExpressionValues(M, theta, eel, tpm, fpkm);
				for (int i = 0; i <= M; i++) {
					params->pme_c[i] += counts[i];
					params->pve_c[i] += double(counts[i]) * counts[i];
					params->pme_tpm[i] += tpm[i];
					params->pme_fpkm[i] += fpkm[i];
				}

				for (int i = 0; i < m; i++) {
				  int b = gi.spAt(i), e = gi.spAt(i + 1);
				  double count = 0.0;
				  for (int j = b; j < e; j++) count += counts[j];
				  params->pve_c_genes[i] += count * count;
				}

				if (alleleS)
				  for (int i = 0; i < m_trans; i++) {
				    int b = ta.spAt(i), e = ta.spAt(i + 1);
				    double count = 0.0;
				    for (int j = b; j < e; j++) count += counts[j];
				    params->pve_c_trans[i] += count * count;
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

	pve_c_genes.assign(m, 0);
	pve_c_trans.clear();
	if (alleleS) pve_c_trans.assign(m_trans, 0);

	for (int i = 0; i < nThreads; i++) {
		fclose(paramsArray[i].fo);
		delete paramsArray[i].engine;
		for (int j = 0; j <= M; j++) {
			pme_c[j] += paramsArray[i].pme_c[j];
			pve_c[j] += paramsArray[i].pve_c[j];
			pme_tpm[j] += paramsArray[i].pme_tpm[j];
			pme_fpkm[j] += paramsArray[i].pme_fpkm[j];
		}

		for (int j = 0; j < m; j++) 
		  pve_c_genes[j] += paramsArray[i].pve_c_genes[j];
		
		if (alleleS) 
		  for (int j = 0; j < m_trans; j++) 
		    pve_c_trans[j] += paramsArray[i].pve_c_trans[j];

		delete[] paramsArray[i].pme_c;
		delete[] paramsArray[i].pve_c;
		delete[] paramsArray[i].pme_tpm;
		delete[] paramsArray[i].pme_fpkm;

		delete[] paramsArray[i].pve_c_genes;
		if (alleleS) delete[] paramsArray[i].pve_c_trans;
	}
	delete[] paramsArray;

	for (int i = 0; i <= M; i++) {
		pme_c[i] /= NSAMPLES;
		pve_c[i] = (pve_c[i] - double(NSAMPLES) * pme_c[i] * pme_c[i]) / double(NSAMPLES - 1);
		if (pve_c[i] < 0.0) pve_c[i] = 0.0;
		pme_tpm[i] /= NSAMPLES;
		pme_fpkm[i] /= NSAMPLES;
	}

	for (int i = 0; i < m; i++) {
	  int b = gi.spAt(i), e = gi.spAt(i + 1);
	  double pme_c_gene = 0.0;
	  for (int j = b; j < e; j++) pme_c_gene += pme_c[j];
	  pve_c_genes[i] = (pve_c_genes[i] - double(NSAMPLES) * pme_c_gene * pme_c_gene) / double(NSAMPLES - 1);
	  if (pve_c_genes[i] < 0.0) pve_c_genes[i] = 0.0;
	}

	if (alleleS) 
	  for (int i = 0; i < m_trans; i++) {
	    int b = ta.spAt(i), e = ta.spAt(i + 1);
	    double pme_c_tran = 0.0;
	    for (int j = b; j < e; j++) pme_c_tran += pme_c[j];
	    pve_c_trans[i] = (pve_c_trans[i] - double(NSAMPLES) * pme_c_tran * pme_c_tran) / double(NSAMPLES - 1);
	    if (pve_c_trans[i] < 0.0) pve_c_trans[i] = 0.0;
	  }
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		// pliu
		// add an option --prior to take priors
		printf("Usage: rsem-run-gibbs reference_name imdName statName BURNIN NSAMPLES GAP [-p #Threads] [--seed seed] [--pseudo-count pseudo_count] [--prior file] [-q]\n");
    printf("\n");
    printf("Format of the prior file:\n");
    printf("- One isoform's prior per line\n");
    printf("- Priors must be in the same order as in the .ti file\n");
    printf("- Priors for those to-be-omitted isoforms must be included as well\n");
    printf("- Comments can be added after prior separated by space(s)\n");
		exit(-1);
	}

	strcpy(refName, argv[1]);
	strcpy(imdName, argv[2]);
	strcpy(statName, argv[3]);

	BURNIN = atoi(argv[4]);
	NSAMPLES = atoi(argv[5]);
	GAP = atoi(argv[6]);

	nThreads = 1;
	hasSeed = false;
	pseudoC = 1.0;
	quiet = false;

	// pliu
	has_prior = false;
	//////

	for (int i = 7; i < argc; i++) {
		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--seed")) {
		  hasSeed = true;
		  int len = strlen(argv[i + 1]);
		  seed = 0;
		  for (int k = 0; k < len; k++) seed = seed * 10 + (argv[i + 1][k] - '0');
		}
		if (!strcmp(argv[i], "--pseudo-count")) pseudoC = atof(argv[i + 1]);
		if (!strcmp(argv[i], "-q")) quiet = true;

		// pliu
		if ( ! strcmp(argv[i], "--prior") ) {
			has_prior = true;
			strcpy(fprior, argv[i+1]);
		}
		//////
	}
	verbose = !quiet;

	assert(NSAMPLES > 1); // Otherwise, we cannot calculate posterior variance

	if (nThreads > NSAMPLES) {
		nThreads = NSAMPLES;
		fprintf(stderr, "Warning: Number of samples is less than number of threads! Change the number of threads to %d!\n", nThreads);
	}

	load_data(refName, statName, imdName);
	load_group_info(refName);
	load_omit_info(imdName);

	// pliu
	// have to do it after load_data() in order to use 'M'
	// the variable 'totc' will be re-calculated by including the prior info
	if ( has_prior ) {
		load_prior_info(fprior);
	}
	//////

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
	
	writeResultsGibbs(M, m, m_trans, gi, gt, ta, alleleS, imdName, pme_c, pme_fpkm, pme_tpm, pve_c, pve_c_genes, pve_c_trans);

	delete mw; // delete the copy

	return 0;
}
