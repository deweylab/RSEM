#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<sstream>
#include<vector>

#include "randomc.h"
#include "utils.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Refs.h"
#include "GroupInfo.h"

using namespace std;

struct Item {
	int sid;
	double conprb;

	Item(int sid, double conprb) {
		this->sid = sid;
		this->conprb = conprb;
	}
};

int model_type;
int m, M, N0, N1, nHits;
double totc;
int BURNIN, CHAINLEN, GAP;
char thetaF[STRLEN], ofgF[STRLEN], groupF[STRLEN], refF[STRLEN], modelF[STRLEN];
char cvsF[STRLEN];

Refs refs;
GroupInfo gi;

vector<double> theta, pme_theta, pme_c, eel;

vector<int> s, z;
vector<Item> hits;
vector<int> counts;

bool quiet;

vector<double> arr;
CRandomMersenne rg(time(NULL));

void load_data(char* reference_name, char* sample_name, char* imdName) {
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
	sprintf(thetaF, "%s.theta",sample_name);
	fin.open(thetaF);
	if (!fin.is_open()) {
		fprintf(stderr, "Cannot open %s!\n", thetaF);
		exit(-1);
	}
	fin>>tmpVal;
	if (tmpVal != M + 1) {
		fprintf(stderr, "Number of transcripts is not consistent in %s and %s!\n", refF, thetaF);
		exit(-1);
	}
	theta.clear(); theta.resize(M + 1);
	for (int i = 0; i <= M; i++) fin>>theta[i];
	fin.close();

	//load ofgF;
	sprintf(ofgF, "%s.ofg", imdName);
	fin.open(ofgF);
	if (!fin.is_open()) {
		fprintf(stderr, "Cannot open %s!\n", ofgF);
		exit(-1);
	}
	fin>>tmpVal>>N0;
	if (tmpVal != M) {
		fprintf(stderr, "M in %s is not consistent with %s!\n", ofgF, refF);
		exit(-1);
	}
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

	if (verbose) { printf("Loading Data is finished!\n"); }
}

// arr should be cumulative!
// interval : [,)
// random number should be in [0, arr[len - 1])
// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0...len-1
int sample(vector<double>& arr, int len) {
  int l, r, mid;
  double prb = rg.Random() * arr[len - 1];

  l = 0; r = len - 1;
  while (l <= r) {
    mid = (l + r) / 2;
    if (arr[mid] <= prb) l = mid + 1;
    else r = mid - 1;
  }

  if (l >= len) { printf("%d %lf %lf\n", len, arr[len - 1], prb); }
  assert(l < len);

  return l;
}

void init() {
	int len, fr, to;

	arr.clear();
	z.clear();
	counts.clear();

	z.resize(N1);
	counts.resize(M + 1, 1); // 1 pseudo count
	counts[0] += N0;

	for (int i = 0; i < N1; i++) {
		fr = s[i]; to = s[i + 1];
		len = to - fr;
		arr.resize(len);
		for (int j = fr; j < to; j++) {
			arr[j - fr] = theta[hits[j].sid] * hits[j].conprb;
			if (j > fr) arr[j - fr] += arr[j - fr - 1];  // cumulative
		}
		z[i] = hits[fr + sample(arr, len)].sid;
		++counts[z[i]];
	}

	totc = N0 + N1 + (M + 1);

	if (verbose) { printf("Initialization is finished!\n"); }
}

void writeCountVector(FILE* fo) {
	for (int i = 0; i < M; i++) {
		fprintf(fo, "%d ", counts[i]);
	}
	fprintf(fo, "%d\n", counts[M]);
}

void Gibbs(char* imdName) {
	FILE *fo;
	int fr, to, len;

	sprintf(cvsF, "%s.countvectors", imdName);
	fo = fopen(cvsF, "w");
	assert(CHAINLEN % GAP == 0);
	fprintf(fo, "%d %d\n", CHAINLEN / GAP, M + 1);
	//fprintf(fo, "%d %d\n", CHAINLEN, M + 1);

	pme_c.clear(); pme_c.resize(M + 1, 0.0);
	pme_theta.clear(); pme_theta.resize(M + 1, 0.0);
	for (int ROUND = 1; ROUND <= BURNIN + CHAINLEN; ROUND++) {

		for (int i = 0; i < N1; i++) {
			--counts[z[i]];
			fr = s[i]; to = s[i + 1]; len = to - fr;
			arr.resize(len);
			for (int j = fr; j < to; j++) {
				arr[j - fr] = counts[hits[j].sid] * hits[j].conprb;
				if (j > fr) arr[j - fr] += arr[j - fr - 1]; //cumulative
			}
			z[i] = hits[fr + sample(arr, len)].sid;
			++counts[z[i]];
		}

		if (ROUND > BURNIN) {
			if ((ROUND - BURNIN -1) % GAP == 0) writeCountVector(fo);
			writeCountVector(fo);
			for (int i = 0; i <= M; i++) { 
			  pme_c[i] += counts[i] - 1;
			  pme_theta[i] += counts[i] / totc;
			}
		}

		if (verbose) { printf("ROUND %d is finished!\n", ROUND); }
	}
	fclose(fo);

	for (int i = 0; i <= M; i++) {
	  pme_c[i] /= CHAINLEN;
	  pme_theta[i] /= CHAINLEN;
	}

	if (verbose) { printf("Gibbs is finished!\n"); }
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
	if (denom <= 0) { fprintf(stderr, "No Expected Effective Length is no less than %.6g?!\n", MINEEL); exit(-1); }
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
	if (denom <= 0) { fprintf(stderr, "No alignable reads?!\n"); exit(-1); }
	//assert(denom > 0);
	for (int i = 1; i <= M; i++) {
		tau[i] /= denom;
	}

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	if (fo == NULL) { fprintf(stderr, "Cannot open %s!\n", outF); exit(-1); }
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", pme_c[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.15g%c", tau[i], (i < M ? '\t' : '\n'));
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	if (fo == NULL) { fprintf(stderr, "Cannot open %s!\n", outF); exit(-1); }
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
		printf("Usage: rsem-run-gibbs reference_name sample_name imdName BURNIN CHAINLEN GAP [-q]\n");
		exit(-1);
	}

	BURNIN = atoi(argv[4]);
	CHAINLEN = atoi(argv[5]);
	GAP = atoi(argv[6]);
	load_data(argv[1], argv[2], argv[3]);

	quiet = false;
	if (argc > 7 && !strcmp(argv[7], "-q")) {
		quiet = true;
	}
	verbose = !quiet;

	init();
	Gibbs(argv[3]);

	sprintf(modelF, "%s.model", argv[2]);
	FILE *fi = fopen(modelF, "r");
	if (fi == NULL) { fprintf(stderr, "Cannot open %s!\n", modelF); exit(-1); }
	fscanf(fi, "%d", &model_type);
	fclose(fi);

	switch(model_type) {
	case 0 : writeEstimatedParameters<SingleModel>(modelF, argv[3]); break;
	case 1 : writeEstimatedParameters<SingleQModel>(modelF, argv[3]); break;
	case 2 : writeEstimatedParameters<PairedEndModel>(modelF, argv[3]); break;
	case 3 : writeEstimatedParameters<PairedEndQModel>(modelF, argv[3]); break;
	}

	return 0;
}
