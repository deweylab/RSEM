#ifndef WRITERESULTS_H_
#define WRITERESULTS_H_

#include<cmath>
#include<cstdio>
#include<vector>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "my_assert.h"
#include "GroupInfo.h"
#include "Transcript.h"
#include "Transcripts.h"
#include "Refs.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

template<class ModelType>
void calcExpectedEffectiveLengths(int M, Refs& refs, ModelType& model, std::vector<double>& eel) {
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
		int pos1 = std::max(std::min(totLen - fullLen + 1, ub) - lb, 0);
		int pos2 = std::max(std::min(totLen, ub) - lb, 0);

		if (pos2 == 0) { eel[i] = 0.0; continue; }
    
		eel[i] = fullLen * cdf[pos1] + ((cdf[pos2] - cdf[pos1]) * (totLen + 1) - (clen[pos2] - clen[pos1]));
		assert(eel[i] >= 0);
		if (eel[i] < MINEEL) { eel[i] = 0.0; }
	}
  
	delete[] pdf;
	delete[] cdf;
	delete[] clen;
}

void polishTheta(int M, std::vector<double>& theta, const std::vector<double>& eel, const double* mw) {
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

void calcExpressionValues(int M, const std::vector<double>& theta, const std::vector<double>& eel, std::vector<double>& tpm, std::vector<double>& fpkm) {
	double denom;
	std::vector<double> frac;

	//calculate fraction of count over all mappabile reads
	denom = 0.0;
	frac.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++) 
	  if (eel[i] >= EPSILON) {
	    frac[i] = theta[i];
	    denom += frac[i];
	  }
	// general_assert(denom >= EPSILON, "No alignable reads?!");
	if (denom < EPSILON) denom = 1.0;
	for (int i = 1; i <= M; i++) frac[i] /= denom;
  
	//calculate FPKM
	fpkm.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++)
		if (eel[i] >= EPSILON) fpkm[i] = frac[i] * 1e9 / eel[i];

	//calculate TPM
	tpm.assign(M + 1, 0.0);
	denom = 0.0;
	for (int i = 1; i <= M; i++) denom += fpkm[i];
	if (denom < EPSILON) denom = 1.0;
	for (int i = 1; i <= M; i++) tpm[i] = fpkm[i] / denom * 1e6;  
}

inline bool isAlleleSpecific(const char* refName, GroupInfo* gt = NULL, GroupInfo* ta = NULL) {
  bool alleleS;
  char gtF[STRLEN], taF[STRLEN];

  sprintf(gtF, "%s.gt", refName);
  sprintf(taF, "%s.ta", refName);
  std::ifstream gtIF(gtF), taIF(taF);
  alleleS = gtIF.is_open() && taIF.is_open();
  if (gtIF.is_open()) gtIF.close();
  if (taIF.is_open()) taIF.close();

  if (alleleS) { 
    if (gt != NULL) gt->load(gtF); 
    if (ta != NULL) ta->load(taF); 
  }

  return alleleS;
}

void writeResultsEM(int M, const char* refName, const char* imdName, Transcripts& transcripts, std::vector<double>& theta, std::vector<double>& eel, double* counts, bool appendNames) {
	char outF[STRLEN];
	FILE *fo;

	int m;
	GroupInfo gi;
	char groupF[STRLEN];

	std::vector<int> tlens;
	std::vector<double> fpkm, tpm, isopct;
	std::vector<double> glens, gene_eels, gene_counts, gene_tpm, gene_fpkm;
	
	// Load group info
	sprintf(groupF, "%s.grp", refName);
	gi.load(groupF);
	m = gi.getm();

	// For allele-specific expression
	int m_trans = 0;
	GroupInfo gt, ta;
	std::vector<double> trans_lens, trans_eels, trans_counts, trans_tpm, trans_fpkm, ta_pct, gt_pct;

	bool alleleS = isAlleleSpecific(refName, &gt, &ta); // if allele-specific

	calcExpressionValues(M, theta, eel, tpm, fpkm);

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

			gene_counts[i] += counts[j];
			gene_tpm[i] += tpm[j];
			gene_fpkm[i] += fpkm[j];
		}

		if (gene_tpm[i] < EPSILON) {
			double frac = 1.0 / (e - b);
			for (int j = b; j < e; j++) {
				glens[i] += tlens[j] * frac;
				gene_eels[i] += eel[j] * frac;
			}
		}
		else {
			for (int j = b; j < e; j++) {
				isopct[j] = gene_tpm[i] > EPSILON ? tpm[j] / gene_tpm[i] : 0.0;
				glens[i] += tlens[j] * isopct[j];
				gene_eels[i] += eel[j] * isopct[j];
			}
		}
	}

	if (alleleS) {
	  m_trans = ta.getm();
	  ta_pct.assign(M + 1, 0.0);
	  trans_lens.assign(m_trans, 0.0); trans_eels.assign(m_trans, 0.0);
	  trans_counts.assign(m_trans, 0.0); trans_tpm.assign(m_trans, 0.0); trans_fpkm.assign(m_trans, 0.0);
	  
	  for (int i = 0; i < m_trans; i++) {
		int b = ta.spAt(i), e = ta.spAt(i + 1);
		for (int j = b; j < e; j++) {
			trans_counts[i] += counts[j];
			trans_tpm[i] += tpm[j];
			trans_fpkm[i] += fpkm[j];
		}

		if (trans_tpm[i] < EPSILON) {
			double frac = 1.0 / (e - b);
			for (int j = b; j < e; j++) {
				trans_lens[i] += tlens[j] * frac;
				trans_eels[i] += eel[j] * frac;
			}
		}
		else {
			for (int j = b; j < e; j++) {
				ta_pct[j] = trans_tpm[i] > EPSILON ? tpm[j] / trans_tpm[i] : 0.0;
				trans_lens[i] += tlens[j] * ta_pct[j];
				trans_eels[i] += eel[j] * ta_pct[j];
			}
		}
	  } 
	  
	  gt_pct.assign(m_trans, 0.0);
	  for (int i = 0; i < m; i++) 
	    if (gene_tpm[i] >= EPSILON) {
	      int b = gt.spAt(i), e = gt.spAt(i + 1);
	      for (int j = b; j < e; j++) gt_pct[j] = gene_tpm[i] > EPSILON ? trans_tpm[j] / gene_tpm[i] : 0.0;
	    }
	}

	if (!alleleS) {
	  //isoform level results
	  sprintf(outF, "%s.iso_res", imdName);
	  fo = fopen(outF, "w");
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);

	    fprintf(fo, "%s", transcript.getTranscriptID().c_str());
	    if (appendNames && transcript.getTranscriptName() != "") 
	      fprintf(fo, "_%s", transcript.getTranscriptName().c_str());
	    fprintf(fo, "%c", (i < M ? '\t' : '\n'));
	  }
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    
	    fprintf(fo, "%s", transcript.getGeneID().c_str());
	    if (appendNames && transcript.getGeneName() != "")
	      fprintf(fo, "_%s", transcript.getGeneName().c_str());
	    fprintf(fo, "%c", (i < M ? '\t' : '\n'));
	  }
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%d%c", tlens[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", eel[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", counts[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", tpm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", fpkm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	  fclose(fo);
	}
	else {
	  // allele level results
	  sprintf(outF, "%s.allele_res", imdName);
	  fo = fopen(outF, "w");
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    fprintf(fo, "%s%c", transcript.getSeqName().c_str(), (i < M ? '\t' : '\n'));
	  }
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    fprintf(fo, "%s%c", transcript.getTranscriptID().c_str(), (i < M ? '\t' : '\n'));
	  }
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    fprintf(fo, "%s%c", transcript.getGeneID().c_str(), (i < M ? '\t' : '\n'));
	  }
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%d%c", tlens[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", eel[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", counts[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", tpm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", fpkm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", ta_pct[i] * 1e2, (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	  fclose(fo);

	  // isoform level results
	  sprintf(outF, "%s.iso_res", imdName);
	  fo = fopen(outF, "w");
	  for (int i = 0; i < m_trans; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(ta.spAt(i));
	    fprintf(fo, "%s%c", transcript.getTranscriptID().c_str(), (i < m_trans - 1 ? '\t' : '\n'));
	  }
	  for (int i = 0; i < m_trans; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(ta.spAt(i));
	    fprintf(fo, "%s%c", transcript.getGeneID().c_str(), (i < m_trans - 1 ? '\t' : '\n'));
	  }
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_lens[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_eels[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_counts[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_tpm[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_fpkm[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", gt_pct[i] * 1e2, (i < m_trans - 1 ? '\t' : '\n'));
	  fclose(fo);
	}

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "w");
	for (int i = 0; i < m; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(gi.spAt(i));
		
		fprintf(fo, "%s", transcript.getGeneID().c_str());
		if (appendNames && transcript.getGeneName() != "")
		  fprintf(fo, "_%s", transcript.getGeneName().c_str());
		fprintf(fo, "%c", (i < m - 1 ? '\t' : '\n'));
	}
	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		std::string curtid = "", tid;
		for (int j = b; j < e; j++) {
			const Transcript& transcript = transcripts.getTranscriptAt(j);
			tid = transcript.getTranscriptID();
			if (curtid != tid) {
			  if (curtid != "") fprintf(fo, ",");
			  fprintf(fo, "%s", tid.c_str());
			  if (appendNames && transcript.getTranscriptName() != "")
			    fprintf(fo, "_%s", transcript.getTranscriptName().c_str());
			  curtid = tid;
			}
		}
		fprintf(fo, "%c", (i < m - 1 ? '\t' : '\n'));
	}
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", glens[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_eels[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_counts[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_tpm[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_fpkm[i], (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	if (verbose) { printf("Expression Results are written!\n"); }
}

void writeResultsGibbs(int M, int m, int m_trans, GroupInfo& gi, GroupInfo &gt, GroupInfo &ta, bool alleleS, char* imdName, std::vector<double>& pme_c, std::vector<double>& pme_fpkm, std::vector<double>& pme_tpm, std::vector<double>& pve_c, std::vector<double>& pve_c_genes, std::vector<double>& pve_c_trans) {
	char outF[STRLEN];
	FILE *fo;

	std::vector<double> isopct;
	std::vector<double> gene_counts, gene_tpm, gene_fpkm;

	// For allele-specific expression
	std::vector<double> trans_counts, trans_tpm, trans_fpkm, ta_pct, gt_pct;

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

	if (alleleS) {
	  ta_pct.assign(M + 1, 0.0);
	  trans_counts.assign(m_trans, 0.0); trans_tpm.assign(m_trans, 0.0); trans_fpkm.assign(m_trans, 0.0);
	  
	  for (int i = 0; i < m_trans; i++) {
		int b = ta.spAt(i), e = ta.spAt(i + 1);
		for (int j = b; j < e; j++) {
			trans_counts[i] += pme_c[j];
			trans_tpm[i] += pme_tpm[j];
			trans_fpkm[i] += pme_fpkm[j];
		}
		if (trans_tpm[i] < EPSILON) continue;
		for (int j = b; j < e; j++) 
		  ta_pct[j] = pme_tpm[j] / trans_tpm[i];	
	  }

	  gt_pct.assign(m_trans, 0.0);
	  for (int i = 0; i < m; i++) 
	    if (gene_tpm[i] >= EPSILON) {
	      int b = gt.spAt(i), e = gt.spAt(i + 1);
	      for (int j = b; j < e; j++) gt_pct[j] = trans_tpm[j] / gene_tpm[i];
	    }
	}

	if (!alleleS) {
	  //isoform level results
	  sprintf(outF, "%s.iso_res", imdName);
	  fo = fopen(outF, "a");
	  general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");
	  
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_c[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++) 
	    fprintf(fo, "%.2f%c", sqrt(pve_c[i]), (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_tpm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_fpkm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	  fclose(fo);	 
	}
	else {
	  //allele level results
	  sprintf(outF, "%s.allele_res", imdName);
	  fo = fopen(outF, "a");
	  general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");
	  
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_c[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++) 
	    fprintf(fo, "%.2f%c", sqrt(pve_c[i]), (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_tpm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", pme_fpkm[i], (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", ta_pct[i] * 1e2, (i < M ? '\t' : '\n'));
	  for (int i = 1; i <= M; i++)
	    fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	  fclose(fo);

	  //isoform level results
	  sprintf(outF, "%s.iso_res", imdName);
	  fo = fopen(outF, "a");
	  general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");
	  
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_counts[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++) 
	    fprintf(fo, "%.2f%c", sqrt(pve_c_trans[i]), (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_tpm[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", trans_fpkm[i], (i < m_trans - 1 ? '\t' : '\n'));
	  for (int i = 0; i < m_trans; i++)
	    fprintf(fo, "%.2f%c", gt_pct[i] * 1e2, (i < m_trans - 1 ? '\t' : '\n'));
	  fclose(fo);
	}
 
	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	general_assert(fo != NULL, "Cannot open " + cstrtos(outF) + "!");

	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.2f%c", gene_counts[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++) 
	  fprintf(fo, "%.2f%c", sqrt(pve_c_genes[i]), (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.2f%c", gene_tpm[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
	  fprintf(fo, "%.2f%c", gene_fpkm[i], (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	if (verbose) { printf("Gibbs based expression values are written!\n"); }
}

void writeResultsSimulation(int M, char* refName, char* outFN, Transcripts& transcripts, std::vector<double>& eel, std::vector<double>& counts) {
	char outF[STRLEN];
	FILE *fo;

	int m;
	GroupInfo gi;
	char groupF[STRLEN];

        // Load group info
        sprintf(groupF, "%s.grp", refName);
        gi.load(groupF);
	m = gi.getm();

	std::vector<int> tlens;
	std::vector<double> tpm, fpkm, isopct;
	std::vector<double> glens, gene_eels, gene_counts, gene_tpm, gene_fpkm;
	
        // For allele-specific expression
	int m_trans = 0;
        GroupInfo gt, ta;
	std::vector<double> trans_lens, trans_eels, trans_counts, trans_tpm, trans_fpkm, ta_pct, gt_pct;

	bool alleleS = isAlleleSpecific(refName, &gt, &ta); // if allele-specific
	
	for (int i = 1; i <= M; i++)
		general_assert(eel[i] > EPSILON || counts[i] <= EPSILON, "An isoform whose effecitve length < " + ftos(MINEEL, 6) + " got sampled!");

	calcExpressionValues(M, counts, eel, tpm, fpkm);

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
	    
	    gene_counts[i] += counts[j];
	    gene_tpm[i] += tpm[j];
	    gene_fpkm[i] += fpkm[j];
	  }

	  if (gene_tpm[i] < EPSILON) {
	    double frac = 1.0 / (e - b);
	    for (int j = b; j < e; j++) {
	      glens[i] += tlens[j] * frac;
	      gene_eels[i] += eel[j] * frac;
	    }
	  }
	  else {
	    for (int j = b; j < e; j++) {
	      isopct[j] = tpm[j] / gene_tpm[i];
	      glens[i] += tlens[j] * isopct[j];
	      gene_eels[i] += eel[j] * isopct[j];
	    }
	  }
        }

        if (alleleS) {
          m_trans = ta.getm();
          ta_pct.assign(M + 1, 0.0);
          trans_lens.assign(m_trans, 0.0); trans_eels.assign(m_trans, 0.0);
          trans_counts.assign(m_trans, 0.0); trans_tpm.assign(m_trans, 0.0); trans_fpkm.assign(m_trans, 0.0);

          for (int i = 0; i < m_trans; i++) {
	    int b = ta.spAt(i), e = ta.spAt(i + 1);
	    for (int j = b; j < e; j++) {
	      trans_counts[i] += counts[j];
	      trans_tpm[i] += tpm[j];
	      trans_fpkm[i] += fpkm[j];
	    }

	    if (trans_tpm[i] < EPSILON) {
	      double frac = 1.0 / (e - b);
	      for (int j = b; j < e; j++) {
		trans_lens[i] += tlens[j] * frac;
		trans_eels[i] += eel[j] * frac;
	      }
	    }
	    else {
	      for (int j = b; j < e; j++) {
		ta_pct[j] = tpm[j] / trans_tpm[i];
		trans_lens[i] += tlens[j] * ta_pct[j];
		trans_eels[i] += eel[j] * ta_pct[j];
	      }
	    }
          }

          gt_pct.assign(m_trans, 0.0);
          for (int i = 0; i < m; i++)
            if (gene_tpm[i] >= EPSILON) {
              int b = gt.spAt(i), e = gt.spAt(i + 1);
              for (int j = b; j < e; j++) gt_pct[j] = trans_tpm[j] / gene_tpm[i];
            }
        }

	//allele level
	if (alleleS) {
	  sprintf(outF, "%s.sim.alleles.results", outFN);
	  fo = fopen(outF, "w");
	  fprintf(fo, "allele_id\ttranscript_id\tgene_id\tlength\teffective_length\tcount\tTPM\tFPKM\tAlleleIsoPct\tAlleleGenePct\n");
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    fprintf(fo, "%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", transcript.getSeqName().c_str(), transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str(), tlens[i],
		    eel[i], counts[i], tpm[i], fpkm[i], ta_pct[i] * 1e2, isopct[i] * 1e2);
	  }
	  fclose(fo);
	}

	//isoform level
	sprintf(outF, "%s.sim.isoforms.results", outFN);
	fo = fopen(outF, "w");
	fprintf(fo, "transcript_id\tgene_id\tlength\teffective_length\tcount\tTPM\tFPKM\tIsoPct\n");
	if (!alleleS) {
	  for (int i = 1; i <= M; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(i);
	    fprintf(fo, "%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str(), tlens[i],
		    eel[i], counts[i], tpm[i], fpkm[i], isopct[i] * 1e2);
	  }
	}
	else {
	  for (int i = 0; i < m_trans; i++) {
	    const Transcript& transcript = transcripts.getTranscriptAt(ta.spAt(i));
	    fprintf(fo, "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str(), trans_lens[i],
                    trans_eels[i], trans_counts[i], trans_tpm[i], trans_fpkm[i], gt_pct[i] * 1e2);
	  }
	}
	fclose(fo);

	//gene level
	sprintf(outF, "%s.sim.genes.results", outFN);
	fo = fopen(outF, "w");
	fprintf(fo, "gene_id\ttranscript_id(s)\tlength\teffective_length\tcount\tTPM\tFPKM\n");
	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		const std::string& gene_id = transcripts.getTranscriptAt(b).getGeneID();
		fprintf(fo, "%s\t", gene_id.c_str());
		std::string curtid = "", tid;
		for (int j = b; j < e; j++) {
			tid = transcripts.getTranscriptAt(j).getTranscriptID();
			if (curtid != tid) {
			  if (curtid != "") fprintf(fo, ",");
			  fprintf(fo, "%s", tid.c_str());
			  curtid = tid;
			}
		}
		fprintf(fo, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", glens[i], gene_eels[i], gene_counts[i], gene_tpm[i], gene_fpkm[i]);
	}
	fclose(fo);
}

#endif
