#ifndef SINGLEMODEL_H_
#define SINGLEMODEL_H_

#include<cmath>
#include<cstdio>
#include<cassert>
#include<cstring>
#include<string>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<vector>

#include "utils.h"
#include "my_assert.h"
#include "Orientation.h"
#include "LenDist.h"
#include "RSPD.h"
#include "Profile.h"
#include "NoiseProfile.h"

#include "ModelParams.h"
#include "RefSeq.h"
#include "Refs.h"
#include "SingleRead.h"
#include "SingleHit.h"
#include "ReadReader.h"

#include "simul.h"

class SingleModel {
public:
	SingleModel(Refs* refs = NULL) {
		this->refs = refs;
		M = (refs != NULL ? refs->getM() : 0);
		memset(N, 0, sizeof(N));
		estRSPD = false;
		needCalcConPrb = true;

		ori = new Orientation();
		gld = new LenDist();
		mld = NULL;
		rspd = new RSPD(estRSPD);
		pro = new Profile();
		npro = new NoiseProfile();

		mean = -1.0; sd = 0.0;
		mw = NULL;

		seedLen = 0;
	}

	//If it is not a master node, only init & update can be used!
	SingleModel(ModelParams& params, bool isMaster = true) {
		M = params.M;
		memcpy(N, params.N, sizeof(params.N));
		refs = params.refs;
		estRSPD = params.estRSPD;
		mean = params.mean; sd = params.sd;
		seedLen = params.seedLen;
		needCalcConPrb = true;

		ori = NULL; gld = NULL; mld = NULL; rspd = NULL; pro = NULL; npro = NULL;
		mw = NULL;

		if (isMaster) {
			gld = new LenDist(params.minL, params.maxL);
			if (mean >= EPSILON) {
				mld = new LenDist(params.mate_minL, params.mate_maxL);
			}
			if (!estRSPD) { rspd = new RSPD(estRSPD); }
		}

		ori = new Orientation(params.probF);
		if (estRSPD) { rspd = new RSPD(estRSPD, params.B); }
		pro = new Profile(params.maxL);
		npro = new NoiseProfile();
	}

	~SingleModel() {
		refs = NULL;
		if (ori != NULL) delete ori;
		if (gld != NULL) delete gld;
		if (mld != NULL) delete mld;
		if (rspd != NULL) delete rspd;
		if (pro != NULL) delete pro;
		if (npro != NULL) delete npro;
		if (mw != NULL) delete[] mw;
		/* delete[] p1, p2 */
	}

	void estimateFromReads(const char*);

	//if prob is too small, just make it 0
	double getConPrb(const SingleRead& read, const SingleHit& hit) {
		if (read.isLowQuality()) return 0.0;

		double prob;
		int sid = hit.getSid();
		RefSeq &ref = refs->getRef(sid);
		int fullLen = ref.getFullLen();
		int totLen = ref.getTotLen();
		int dir = hit.getDir();
		int pos = hit.getPos();
		int readLen = read.getReadLength();
		int fpos = (dir == 0 ? pos : totLen - pos - readLen); // the aligned position reported in SAM file, should be a coordinate in forward strand

		general_assert(fpos >= 0, "The alignment of read " + read.getName() + " to transcript " + itos(sid) + " starts at " + itos(fpos) + \
				" from the forward direction, which should be a non-negative number! " + \
				"It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(fpos + readLen <= totLen,"Read " + read.getName() + " is hung over the end of transcript " + itos(sid) + "! " \
				+ "It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(readLen <= totLen, "Read " + read.getName() + " has length " + itos(readLen) + ", but it is aligned to transcript " \
				+ itos(sid) + ", whose length (" + itos(totLen) + ") is shorter than the read's length!");

		int seedPos = (dir == 0 ? pos : totLen - pos - seedLen); // the aligned position of the seed in forward strand coordinates
		if (seedPos >= fullLen || ref.getMask(seedPos)) return 0.0;

		int effL;
		double value;

		if (mld != NULL) {
			int minL = std::max(readLen, gld->getMinL());
			int maxL = std::min(totLen - pos, gld->getMaxL());
			int pfpos; // possible fpos for fragment
			value = 0.0;
			for (int fragLen = minL; fragLen <= maxL; fragLen++) {
				pfpos = (dir == 0 ? pos : totLen - pos - fragLen);
				effL = std::min(fullLen, totLen - fragLen + 1);
				value += gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * mld->getAdjustedProb(readLen, fragLen);
			}
		}
		else {
			effL = std::min(fullLen, totLen - readLen + 1);
			value = gld->getAdjustedProb(readLen, totLen) * rspd->getAdjustedProb(fpos, effL, fullLen);
		}

		prob = ori->getProb(dir) * value * pro->getProb(read.getReadSeq(), ref, pos, dir);

		if (prob < EPSILON) { prob = 0.0; }


		prob = (mw[sid] < EPSILON ? 0.0 : prob / mw[sid]);

		return prob;
	}

	double getNoiseConPrb(const SingleRead& read) {
		if (read.isLowQuality()) return 0.0;
		double prob = mld != NULL ? mld->getProb(read.getReadLength()) : gld->getProb(read.getReadLength());
		prob *= npro->getProb(read.getReadSeq());
		if (prob < EPSILON) { prob = 0.0; }

		prob = (mw[0] < EPSILON ? 0.0 : prob / mw[0]);

		return prob;
	}

	double getLogP() { return npro->getLogP(); }

	void init();

	void update(const SingleRead& read, const SingleHit& hit, double frac) {
		if (read.isLowQuality() || frac < EPSILON) return;

		RefSeq& ref = refs->getRef(hit.getSid());
		int dir = hit.getDir();
		int pos = hit.getPos();

		if (estRSPD) {
			int fullLen = ref.getFullLen();

			// Only use one strand to estimate RSPD
			if (ori->getProb(0) >= ORIVALVE && dir == 0) {
				rspd->update(pos, fullLen, frac);
			}

			if (ori->getProb(0) < ORIVALVE && dir == 1) {
				int totLen = ref.getTotLen();
				int readLen = read.getReadLength();

				int pfpos, effL; 

				if (mld != NULL) {
					int minL = std::max(readLen, gld->getMinL());
					int maxL = std::min(totLen - pos, gld->getMaxL());
					double sum = 0.0;
					assert(maxL >= minL);
					std::vector<double> frag_vec(maxL - minL + 1, 0.0);

					for (int fragLen = minL; fragLen <= maxL; fragLen++) {
						pfpos = totLen - pos - fragLen;
						effL = std::min(fullLen, totLen - fragLen + 1);
						frag_vec[fragLen - minL] = gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * mld->getAdjustedProb(readLen, fragLen);
						sum += frag_vec[fragLen - minL];
					}
					assert(sum >= EPSILON);
					for (int fragLen = minL; fragLen <= maxL; fragLen++) {
						pfpos = totLen - pos - fragLen;
						rspd->update(pfpos, fullLen, frac * (frag_vec[fragLen - minL] / sum));
					}
				}
				else {
					rspd->update(totLen - pos - readLen, fullLen, frac);
				}
			}
		}
		pro->update(read.getReadSeq(), ref, pos, dir, frac);
	}

	void updateNoise(const SingleRead& read, double frac) {
		if (read.isLowQuality() || frac < EPSILON) return;

		npro->update(read.getReadSeq(), frac);
	}

	void finish();

	void collect(const SingleModel&);

	bool getNeedCalcConPrb() { return needCalcConPrb; }
	void setNeedCalcConPrb(bool value) { needCalcConPrb = value; }

	//void calcP1();
	//void calcP2();
	//double* getP1() { return p1; }
	//double* getP2() { return p2; }

	void read(const char*);
	void write(const char*);

	const LenDist& getGLD() { return *gld; }

	void startSimulation(simul*, const std::vector<double>&);
	bool simulate(READ_INT_TYPE, SingleRead&, int&);
	void finishSimulation();

	const double* getMW() { 
	  assert(mw != NULL);
	  return mw;
	}

	int getModelType() const { return model_type; }

private:
	static const int model_type = 0;
	static const int read_type = 0;

	int M;
	READ_INT_TYPE N[3];
	Refs *refs;
	double mean, sd;
	int seedLen;
	//double *p1, *p2; P_i' & P_i''

	bool estRSPD; // true if estimate RSPD
	bool needCalcConPrb; // true need, false does not need

	Orientation *ori;
	LenDist *gld, *mld;
	RSPD *rspd;
	Profile *pro;
	NoiseProfile *npro;

	simul *sampler; // for simulation
	double *theta_cdf; // for simulation

	double *mw; // for masking

	void calcMW();
};

void SingleModel::estimateFromReads(const char* readFN) {
	int s;
	char readFs[2][STRLEN];
	SingleRead read;

	int n_warns = 0;
	
	mld != NULL ? mld->init() : gld->init();
	
	for (int i = 0; i < 3; i++)
		if (N[i] > 0) {
			genReadFileNames(readFN, i, read_type, s, readFs);
			ReadReader<SingleRead> reader(s, readFs, refs->hasPolyA(), seedLen); // allow calculation of calc_lq() function

			READ_INT_TYPE cnt = 0;
			while (reader.next(read)) {
				if (!read.isLowQuality()) {
					mld != NULL ? mld->update(read.getReadLength(), 1.0) : gld->update(read.getReadLength(), 1.0);
					if (i == 0) { npro->updateC(read.getReadSeq()); }
				}
				else if (read.getReadLength() < seedLen)
				  if (++n_warns <= MAX_WARNS)
				    fprintf(stderr, "Warning: Read %s is ignored due to read length (= %d) < seed length (= %d)!\n", read.getName().c_str(), read.getReadLength(), seedLen);
								
				++cnt;
				if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " READS PROCESSED"<< std::endl; }
			}

			if (verbose) { std::cout<< "estimateFromReads, N"<< i<< " finished."<< std::endl; }
		}

	if (n_warns > 0) fprintf(stderr, "Warning: There are %d reads ignored in total.\n", n_warns);
	
	mld != NULL ? mld->finish() : gld->finish();
	if (mean >= EPSILON) { //mean should be > 0
	  assert(mld->getMaxL() <= gld->getMaxL());
	  gld->setAsNormal(mean, sd, std::max(mld->getMinL(), gld->getMinL()), gld->getMaxL());
	}
	npro->calcInitParams();

	mw = new double[M + 1];
	calcMW();
}

void SingleModel::init() {
	if (estRSPD) rspd->init();
	pro->init();
	npro->init();
}

void SingleModel::finish() {
	if (estRSPD) rspd->finish();
	pro->finish();
	npro->finish();
	needCalcConPrb = true;
	if (estRSPD) calcMW();
}

void SingleModel::collect(const SingleModel& o) {
	if (estRSPD) rspd->collect(*(o.rspd));
	pro->collect(*(o.pro));
	npro->collect(*(o.npro));
}

//Only master node can call
void SingleModel::read(const char* inpF) {
	int val;
	FILE *fi = fopen(inpF, "r");

	general_assert(fi != NULL, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	assert(fscanf(fi, "%d", &val) == 1);
	assert(val == model_type);

	ori->read(fi);
	gld->read(fi);
	assert(fscanf(fi, "%d", &val) == 1);
	if (val > 0) {
		if (mld == NULL) mld = new LenDist();
		mld->read(fi);
	}
	rspd->read(fi);
	pro->read(fi);
	npro->read(fi);

	if (fscanf(fi, "%d", &val) == 1) {
		if (M == 0) M = val;
		if (M == val) {
			mw = new double[M + 1];
			for (int i = 0; i <= M; i++) assert(fscanf(fi, "%lf", &mw[i]) == 1);
		}
	}

	fclose(fi);
}

//Only master node can call. Only be called at EM.cpp
void SingleModel::write(const char* outF) {
	FILE *fo = fopen(outF, "w");

	fprintf(fo, "%d\n", model_type);
	fprintf(fo, "\n");

	ori->write(fo);  fprintf(fo, "\n");
	gld->write(fo);  fprintf(fo, "\n");
	if (mld != NULL) {
		fprintf(fo, "1\n");
		mld->write(fo);
	}
	else { fprintf(fo, "0\n"); }
	fprintf(fo, "\n");
	rspd->write(fo); fprintf(fo, "\n");
	pro->write(fo);  fprintf(fo, "\n");
	npro->write(fo);

	if (mw != NULL) {
	  fprintf(fo, "\n%d\n", M);
	  for (int i = 0; i < M; i++) {
	    fprintf(fo, "%.15g ", mw[i]);
	  }
	  fprintf(fo, "%.15g\n", mw[M]);
	}

	fclose(fo);
}

void SingleModel::startSimulation(simul* sampler, const std::vector<double>& theta) {
	this->sampler = sampler;

	theta_cdf = new double[M + 1];
	for (int i = 0; i <= M; i++) {
		theta_cdf[i] = theta[i];
		if (i > 0) theta_cdf[i] += theta_cdf[i - 1];
	}

	rspd->startSimulation(M, refs);
	pro->startSimulation();
	npro->startSimulation();
}

bool SingleModel::simulate(READ_INT_TYPE rid, SingleRead& read, int& sid) {
	int dir, pos, readLen, fragLen;
	std::string name;
	std::string readseq;
	std::ostringstream strout;

	sid = sampler->sample(theta_cdf, M + 1);

	if (sid == 0) {
		dir = pos = 0;
		readLen = (mld != NULL ? mld->simulate(sampler, -1) : gld->simulate(sampler, -1));
		readseq = npro->simulate(sampler, readLen);
	}
	else {
		RefSeq &ref = refs->getRef(sid);
		dir = ori->simulate(sampler);
		fragLen = gld->simulate(sampler, ref.getTotLen());
		if (fragLen < 0) return false;
		int effL = std::min(ref.getFullLen(), ref.getTotLen() - fragLen + 1);
		pos = rspd->simulate(sampler, sid, effL);
		if (pos < 0) return false;
		if (dir > 0) pos = ref.getTotLen() - pos - fragLen;

		if (mld != NULL) {
			readLen = mld->simulate(sampler, fragLen);
			if (readLen < 0) return false;
			readseq = pro->simulate(sampler, readLen, pos, dir, ref);
		}
		else {
			readseq = pro->simulate(sampler, fragLen, pos, dir, ref);
		}
	}

	strout<<rid<<"_"<<dir<<"_"<<sid<<"_"<<pos;
	name = strout.str();

	read = SingleRead(name, readseq);

	return true;
}

void SingleModel::finishSimulation() {
	delete[] theta_cdf;

	rspd->finishSimulation();
	pro->finishSimulation();
	npro->finishSimulation();
}

void SingleModel::calcMW() {
	double probF, probR;

	assert((mld == NULL ? gld->getMinL() : mld->getMinL()) >= seedLen);
  
	memset(mw, 0, sizeof(double) * (M + 1));
	mw[0] = 1.0;

	probF = ori->getProb(0);
	probR = ori->getProb(1);

	for (int i = 1; i <= M; i++) {
		RefSeq& ref = refs->getRef(i);
		int totLen = ref.getTotLen();
		int fullLen = ref.getFullLen();
		double value = 0.0;
		int minL, maxL;
		int effL, pfpos;
		int end = std::min(fullLen, totLen - seedLen + 1);
		double factor;

		for (int seedPos = 0; seedPos < end; seedPos++)
			if (ref.getMask(seedPos)) {
				//forward
				minL = gld->getMinL();
				maxL = std::min(gld->getMaxL(), totLen - seedPos);
				pfpos = seedPos;
				for (int fragLen = minL; fragLen <= maxL; fragLen++) {
					effL = std::min(fullLen, totLen - fragLen + 1);
					factor = (mld == NULL ? 1.0 : mld->getAdjustedCumulativeProb(std::min(mld->getMaxL(), fragLen), fragLen));
					value += probF * gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * factor;
				}
				//reverse
				minL = gld->getMinL();
				maxL = std::min(gld->getMaxL(), seedPos + seedLen);
				for (int fragLen = minL; fragLen <= maxL; fragLen++) {
					pfpos = seedPos - (fragLen - seedLen);
					effL = std::min(fullLen, totLen - fragLen + 1);
					factor = (mld == NULL ? 1.0 : mld->getAdjustedCumulativeProb(std::min(mld->getMaxL(), fragLen), fragLen));
					value += probR * gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * factor;
				}
			}
    
		//for reverse strand masking
		for (int seedPos = end; seedPos <= totLen - seedLen; seedPos++) {
			minL = std::max(gld->getMinL(), seedPos + seedLen - fullLen + 1);
			maxL = std::min(gld->getMaxL(), seedPos + seedLen);
			for (int fragLen = minL; fragLen <= maxL; fragLen++) {
				pfpos = seedPos - (fragLen - seedLen);
				effL = std::min(fullLen, totLen - fragLen + 1);
				factor = (mld == NULL ? 1.0 : mld->getAdjustedCumulativeProb(std::min(mld->getMaxL(), fragLen), fragLen));
				value += probR * gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * factor;
			}
		}

		mw[i] = 1.0 - value;

		if (mw[i] < 1e-8) {
			//      fprintf(stderr, "Warning: %dth reference sequence is masked for almost all positions!\n", i);
			mw[i] = 0.0;
		}
	}
}

#endif /* SINGLEMODEL_H_ */
