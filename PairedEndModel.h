#ifndef PAIREDENDMODEL_H_
#define PAIREDENDMODEL_H_

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
#include "PairedEndRead.h"
#include "PairedEndHit.h"
#include "ReadReader.h"

#include "simul.h"

class PairedEndModel {
public:
	PairedEndModel(Refs* refs = NULL) {
		this->refs = refs;
		M = (refs != NULL ? refs->getM() : 0);
		memset(N, 0, sizeof(N));
		estRSPD = false;
		needCalcConPrb = true;

		ori = new Orientation();
		gld = new LenDist();
		rspd = new RSPD(estRSPD);
		pro = new Profile();
		npro = new NoiseProfile();
		mld = new LenDist();

		mw = NULL;
		seedLen = 0;
	}

	//If it is not a master node, only init & update can be used!
	PairedEndModel(ModelParams& params, bool isMaster = true) {
		M = params.M;
		memcpy(N, params.N, sizeof(params.N));
		refs = params.refs;
		estRSPD = params.estRSPD;
		seedLen = params.seedLen;
		needCalcConPrb = true;

		ori = NULL; gld = NULL; rspd = NULL; pro = NULL; npro = NULL; mld = NULL;
		mw = NULL;

		if (isMaster) {
			if (!estRSPD) rspd = new RSPD(estRSPD);
			mld = new LenDist(params.mate_minL, params.mate_maxL);
		}

		ori = new Orientation(params.probF);
		gld = new LenDist(params.minL, params.maxL);
		if (estRSPD) rspd = new RSPD(estRSPD, params.B);
		pro = new Profile(params.maxL);
		npro = new NoiseProfile();
	}

	~PairedEndModel() {
		refs = NULL;
		if (ori != NULL) delete ori;
		if (gld != NULL) delete gld;
		if (rspd != NULL) delete rspd;
		if (pro != NULL) delete pro;
		if (npro != NULL) delete npro;
		if (mld != NULL) delete mld;
		if (mw != NULL) delete mw;
	}

	void estimateFromReads(const char*);

	//if prob is too small, just make it 0
	double getConPrb(const PairedEndRead& read, const PairedEndHit& hit) {
		if (read.isLowQuality()) return 0.0;

		double prob;
		int sid = hit.getSid();
		RefSeq &ref = refs->getRef(sid);
		int dir = hit.getDir();
		int pos = hit.getPos();
		int fullLen = ref.getFullLen();
		int totLen = ref.getTotLen();
		int insertLen = hit.getInsertL();

		int fpos = (dir == 0 ? pos : totLen - pos - insertLen); // the aligned position reported in SAM file, should be a coordinate in forward strand
		int effL = std::min(fullLen, totLen - insertLen + 1);
		
		general_assert(fpos >= 0, "The alignment of fragment " + read.getName() + " to transcript " + itos(sid) + " starts at " + itos(fpos) + \
				" from the forward direction, which should be a non-negative number! " + \
				"It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(fpos + insertLen <= totLen,"Fragment " + read.getName() + " is hung over the end of transcript " + itos(sid) + "! " \
				+ "It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(insertLen <= totLen, "Fragment " + read.getName() + " has length " + itos(insertLen) + ", but it is aligned to transcript " \
				+ itos(sid) + ", whose length (" + itos(totLen) + ") is shorter than the fragment's length!");


		if (fpos >= fullLen || ref.getMask(fpos)) return 0.0; // For paired-end model, fpos is the seedPos

		prob = ori->getProb(dir) * gld->getAdjustedProb(insertLen, totLen) *
		       rspd->getAdjustedProb(fpos, effL, fullLen);

		const SingleRead& mate1 = read.getMate1();
		prob *= mld->getAdjustedProb(mate1.getReadLength(), insertLen) *
		        pro->getProb(mate1.getReadSeq(), ref, pos, dir);

		const SingleRead& mate2 = read.getMate2();
		int m2pos = totLen - pos - insertLen;
		int m2dir = !dir;
		prob *= mld->getAdjustedProb(mate2.getReadLength(), insertLen) *
		        pro->getProb(mate2.getReadSeq(), ref, m2pos, m2dir);

		if (prob < EPSILON) { prob = 0.0; }

		prob = (mw[sid] < EPSILON ? 0.0 : prob / mw[sid]);

		return prob;
	}

	double getNoiseConPrb(const PairedEndRead& read) {
		if (read.isLowQuality()) return 0.0;
		double prob;
		const SingleRead& mate1 = read.getMate1();
		const SingleRead& mate2 = read.getMate2();

		prob = mld->getProb(mate1.getReadLength()) * npro->getProb(mate1.getReadSeq());
		prob *= mld->getProb(mate2.getReadLength()) * npro->getProb(mate2.getReadSeq());

		if (prob < EPSILON) { prob = 0.0; }

		prob = (mw[0] < EPSILON ? 0.0: prob / mw[0]);

		return prob;
	}

	double getLogP() { return npro->getLogP(); }

	void init();

	void update(const PairedEndRead& read, const PairedEndHit& hit, double frac) {
		if (read.isLowQuality() || frac < EPSILON) return;

		RefSeq& ref = refs->getRef(hit.getSid());
		const SingleRead& mate1 = read.getMate1();
		const SingleRead& mate2 = read.getMate2();

		gld->update(hit.getInsertL(), frac);
		if (estRSPD) {
			int fpos = (hit.getDir() == 0 ? hit.getPos() : ref.getTotLen() - hit.getPos() - hit.getInsertL());
			rspd->update(fpos, ref.getFullLen(), frac);
		}
		pro->update(mate1.getReadSeq(), ref, hit.getPos(), hit.getDir(), frac);

		int m2pos = ref.getTotLen() - hit.getPos() - hit.getInsertL();
		int m2dir = !hit.getDir();
		pro->update(mate2.getReadSeq(), ref, m2pos, m2dir, frac);
	}

	void updateNoise(const PairedEndRead& read, double frac) {
		if (read.isLowQuality() || frac < EPSILON) return;

		const SingleRead& mate1 = read.getMate1();
		const SingleRead& mate2 = read.getMate2();

		npro->update(mate1.getReadSeq(), frac);
		npro->update(mate2.getReadSeq(), frac);
	}

	void finish();

	void collect(const PairedEndModel&);

	bool getNeedCalcConPrb() { return needCalcConPrb; }
	void setNeedCalcConPrb(bool value) { needCalcConPrb = value; }

	void read(const char*);
	void write(const char*);

	const LenDist& getGLD() { return *gld; }

	void startSimulation(simul*, const std::vector<double>&);
	bool simulate(READ_INT_TYPE, PairedEndRead&, int&);
	void finishSimulation();

	//Use it after function 'read' or 'estimateFromReads'
	const double* getMW() { 
	  assert(mw != NULL);
	  return mw;
	}

	int getModelType() const { return model_type; }

private:
	static const int model_type = 2;
	static const int read_type = 2;

	int M;
	READ_INT_TYPE N[3];
	Refs *refs;
	int seedLen;

	bool estRSPD;
	bool needCalcConPrb; //true need, false does not need

	Orientation *ori;
	LenDist *gld, *mld; //mld1 mate_length_dist
	RSPD *rspd;
	Profile *pro;
	NoiseProfile *npro;

	simul *sampler; // for simulation
	double *theta_cdf; // for simulation

	double *mw; // for masking

	void calcMW();
};

void PairedEndModel::estimateFromReads(const char* readFN) {
    int s;
    char readFs[2][STRLEN];
    PairedEndRead read;

    int n_warns = 0;
    
    mld->init();
    for (int i = 0; i < 3; i++)
    	if (N[i] > 0) {
    		genReadFileNames(readFN, i, read_type, s, readFs);
    		ReadReader<PairedEndRead> reader(s, readFs, refs->hasPolyA(), seedLen); // allow calculation of calc_lq() function

    		READ_INT_TYPE cnt = 0;
    		while (reader.next(read)) {
    			SingleRead mate1 = read.getMate1();
    			SingleRead mate2 = read.getMate2();
			
    			if (!read.isLowQuality()) {
    				mld->update(mate1.getReadLength(), 1.0);
    				mld->update(mate2.getReadLength(), 1.0);
			  
    				if (i == 0) {
    					npro->updateC(mate1.getReadSeq());
    					npro->updateC(mate2.getReadSeq());
    				}
    			}
    			else if (mate1.getReadLength() < seedLen || mate2.getReadLength() < seedLen)
			  if (++n_warns <= MAX_WARNS)
			    fprintf(stderr, "Warning: Read %s is ignored due to at least one of the mates' length < seed length (= %d)!\n", read.getName().c_str(), seedLen);

    			++cnt;
    			if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " READS PROCESSED"<< std::endl; }
    		}

    		if (verbose) { std::cout<< "estimateFromReads, N"<< i<< " finished."<< std::endl; }
    	}

    if (n_warns > 0) fprintf(stderr, "Warning: There are %d reads ignored in total.\n", n_warns);
    
    mld->finish();
    npro->calcInitParams();

    mw = new double[M + 1];
    calcMW();
}

void PairedEndModel::init() {
	gld->init();
	if (estRSPD) rspd->init();
	pro->init();
	npro->init();
}

void PairedEndModel::finish() {
	gld->finish();
	if (estRSPD) rspd->finish();
	pro->finish();
	npro->finish();
	needCalcConPrb = true;
	calcMW();
}

void PairedEndModel::collect(const PairedEndModel& o) {
	gld->collect(*(o.gld));
	if (estRSPD) rspd->collect(*(o.rspd));
	pro->collect(*(o.pro));
	npro->collect(*(o.npro));
}

//Only master node can call
void PairedEndModel::read(const char* inpF) {
	int val;
	FILE *fi = fopen(inpF, "r");

	general_assert(fi != NULL, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	assert(fscanf(fi, "%d", &val) == 1);
	assert(val == model_type);

	ori->read(fi);
	gld->read(fi);
	mld->read(fi);
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
void PairedEndModel::write(const char* outF) {
	FILE *fo = fopen(outF, "w");

	fprintf(fo, "%d\n", model_type);
	fprintf(fo, "\n");

	ori->write(fo);  fprintf(fo, "\n");
	gld->write(fo);  fprintf(fo, "\n");
	mld->write(fo);  fprintf(fo, "\n");
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

void PairedEndModel::startSimulation(simul* sampler, const std::vector<double>& theta) {
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

bool PairedEndModel::simulate(READ_INT_TYPE rid, PairedEndRead& read, int& sid) {
	int dir, pos;
	int insertL, mateL1, mateL2;
	std::string name;
	std::string readseq1, readseq2;
	std::ostringstream strout;

	sid = sampler->sample(theta_cdf, M + 1);

	if (sid == 0) {
		dir = pos = insertL = 0;
		mateL1 = mld->simulate(sampler, -1);
		readseq1 = npro->simulate(sampler, mateL1);

		mateL2 = mld->simulate(sampler, -1);
		readseq2 = npro->simulate(sampler, mateL2);
	}
	else {
		RefSeq &ref = refs->getRef(sid);
		dir = ori->simulate(sampler);
		insertL = gld->simulate(sampler, ref.getTotLen());
		if (insertL < 0) return false;
		int effL = std::min(ref.getFullLen(), ref.getTotLen() - insertL + 1);
		pos = rspd->simulate(sampler, sid, effL);
		if (pos < 0) return false;
		if (dir > 0) pos = ref.getTotLen() - pos - insertL;

		mateL1 = mld->simulate(sampler, insertL);
		readseq1 = pro->simulate(sampler, mateL1, pos, dir, ref);

		int m2pos = ref.getTotLen() - pos - insertL;
		int m2dir = !dir;

		mateL2 = mld->simulate(sampler, insertL);
		readseq2 = pro->simulate(sampler, mateL2, m2pos, m2dir, ref);
	}

	strout<<rid<<"_"<<dir<<"_"<<sid<<"_"<<pos<<"_"<<insertL;
	name = strout.str();

	read = PairedEndRead(SingleRead(name + "/1", readseq1), SingleRead(name + "/2", readseq2));

	return true;
}

void PairedEndModel::finishSimulation() {
	delete[] theta_cdf;

	rspd->finishSimulation();
	pro->finishSimulation();
	npro->finishSimulation();
}

void PairedEndModel::calcMW() {
	assert(mld->getMinL() >= seedLen);

	memset(mw, 0, sizeof(double) * (M + 1));
	mw[0] = 1.0;

	for (int i = 1; i <= M; i++) {
		RefSeq& ref = refs->getRef(i);
		int totLen = ref.getTotLen();
		int fullLen = ref.getFullLen();
		int end = std::min(fullLen, totLen - gld->getMinL() + 1);
		double value = 0.0;
		int minL, maxL;
		int effL, pfpos;

		//seedPos is fpos here
		for (int seedPos = 0; seedPos < end; seedPos++)
			if (ref.getMask(seedPos)) {
				minL = gld->getMinL();
				maxL = std::min(gld->getMaxL(), totLen - seedPos);
				pfpos = seedPos;
				for (int fragLen = minL; fragLen <= maxL; fragLen++) {
					effL = std::min(fullLen, totLen - fragLen + 1);
					value += gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen);
				}
			}

		mw[i] = 1.0 - value;

		if (mw[i] < 1e-8) {
			//fprintf(stderr, "Warning: %dth reference sequence is masked for almost all positions!\n", i);
			mw[i] = 0.0;
		}
	}
}

#endif /* PAIREDENDMODEL_H_ */
