/* Copyright (c) 2017
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#ifndef FRAGLENDIST_H_
#define FRAGLENDIST_H_

#include <vector>

class FragLenDist {
public:
	FragLenDist(int minL, int maxL);

	void setUpCounts();
	void setUpTrans(int M, int *lens, double *efflens, double *countTrans);
	
	int getMinL() const { return lb + 1; }
	int getMaxL() const { return ub; }

	// (len - ub) * cdf[span] because when we perform projection gradient, cdf[span] is not necessarily 1.
	double getEffLen(int len) const { return (len <= lb ? 0 : (len <= ub ? psum[len - lb] : psum[span] + (len - ub) * cdf[span])); }

	void optimize();

	bool operator() (int a, int b) const {
		return pmf[a] < pmf[b];
	}
	
private:
	int lb, ub, span; // (lb, ub], span = ub - lb
	std::vector<double> pmf, cdf, psum; // pmf, probability mass function; cdf, cumulative density function; psum[i] = \sum_{j=1}^{i} (i - j + 1) pmf[j]

	static const double max_step_size, alpha, beta;
	double t;
	std::vector<double> pmf_orig; // temporary vector used to store the original pmf vector during back tracking
	
	int l_lower_SE; // the mininum fragment length - 1 - lb allowed for single-end reads
	bool hasSE, hasPE;
	std::vector<double> countUpper; // for single-end reads, count the number of reads by their maximum fragment lengths; 
	std::vector<double> countFrags; // for paired-end reads, count the number of reads falling into each fragment length bin

	std::vector<double> gd, gd_buffer; // gd, gradient; gd_buffer, buffer for calculating gradient with the efflens
	std::vector<int> ids;
	
	int M;
	int* lens; // each transcript's length
	double* efflens; // each transcript's effective length
	double* countTrans; // number of reads falling into each transcript
	
	void updateAux();

	/*
		@return   the log likelihood function
	*/
	double calc_f();

	/*
		@return   ||the gradient||^2_2
	*/
	double calc_g(std::vector<double>& gd);

	void backtracking();
	
	void projection();
	
};



#endif
