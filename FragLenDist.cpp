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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>

#include "FragLenDist.hpp"

const double FragLenDist::max_step_size = 1.0;
const double FragLenDist::alpha = 0.3;
const double FragLenDist::beta = 0.5;

FragLenDist::FragLenDist(int minL, int maxL) {
	lb = minL - 1;
	ub = maxL;
	span = ub - lb;

	hasSE = hasPE = false;
	countUpper.assign(span + 1, 0.0);
	countFrags.assign(span + 1, 0.0);

	M = 0; lens = NULL; efflens = NULL; countTrans = NULL;
}

void FragLenDist::setUpTrans(int M, int *lens, double *efflens, double *countTrans) {
	this->M = M;
	this->lens = lens;
	this->efflens = efflens;
	this->countTrans = countTrans;

	pmf.assign(span + 1, 0.0);
	cdf.assign(span + 1, 0.0);
	psum.assign(span + 1, 0.0);
}

void FragLenDist::setUpCounts() {
	hasSE = true;
	l_lower_SE = 0;
	countUpper[1] = 1;
	countUpper[2] = 5;
	countUpper[3] = 23;
	//  hasPE = true;
	//  countFrags[1] = 9;
	//  countFrags[2] = 14.5;
	//  countFrags[3] = 5.5;
}

void FragLenDist::optimize() {
	t = 1.0; // initial step size is 1.0
	
	ids.assign(span, 0);
	
	// Uniform dist
	for (int i = 1; i <= span; ++i) pmf[i] = 1.0 / span;
	updateAux();
	printf("f value = %.10g\n", calc_f());
	
	for (int i = 0; i < 10; ++i) {
		backtracking();
		projection();
		updateAux();
		for (int j = 1; j <= span; ++j) printf("%.10g\t", pmf[j]); printf("\n");
		printf("f value = %.10g\n", calc_f());
	}
}

void FragLenDist::updateAux() {
	cdf[0] = psum[0] = 0.0;
	for (int i = 1; i <= span; ++i) {
		cdf[i] = cdf[i - 1] + pmf[i];
		psum[i] = psum[i - 1] + cdf[i];
	}
	for (int i = 1; i <= M; ++i) efflens[i] = getEffLen(lens[i]);
}

double FragLenDist::calc_f() {
	double res = 0.0;
	
	if (hasSE)
		for (int i = l_lower_SE + 1; i <= span; ++i)
			if (countUpper[i] > 0.0) res += countUpper[i] * log(cdf[i] - cdf[l_lower_SE]);
	
	if (hasPE)
		for (int i = 1; i <= span; ++i)
			if (countFrags[i] > 0.0) res += countFrags[i] * log(pmf[i]);
	
	for (int i = 1; i <= M; ++i)
		if (countTrans[i] > 0.0) res -= countTrans[i] * log(efflens[i]);
	
	return res;
}

double FragLenDist::calc_g(std::vector<double>& gd) {
	int id;
	double value;
	gd.assign(span + 1, 0.0); // gradient starts from 1
	
	gd_buffer.assign(span + 1, 0.0);
	for (int i = 1; i <= M; ++i)
		if (countTrans[i] > 0.0) {
			value = countTrans[i] / efflens[i];
			if (lens[i] > ub) {
	gd[span] -= value * (lens[i] - ub + 1);
	gd_buffer[span] -= value;
			}
			else if (lens[i] > lb) {
	id = lens[i] - lb;
	gd[id] -= value;
	gd_buffer[id] -= value;
			}
		}
	
	for (int i = span - 1; i > 0; --i) {
		gd[i] += gd[i + 1] + gd_buffer[i + 1];
		gd_buffer[i] += gd_buffer[i + 1];
	}
	
	if (hasSE) {
		double delta = 0.0;
		for (int i = span; i > l_lower_SE; --i)
			if (countUpper[i] > 0.0) {
	delta += countUpper[i] / (cdf[i] - cdf[l_lower_SE]);
	gd[i] += delta;
			}
	}
	
	if (hasPE)
		for (int i = 1; i <= span; ++i)
			if (countFrags[i] > 0.0) gd[i] += countFrags[i] / pmf[i];

	value = 0.0;
	for (int i = 1; i <= span; ++i) value += gd[i] * gd[i];
	return value;
}

// Implemented a not so naive back tracking scheme, the resulting pmf vector might need call updateAux()
void FragLenDist::backtracking() {
	int status; // 0, first try; 1, t *= beta; 2, t /= beta
	bool noNeg, success; // noNeg, no negative values
	double fvalue, new_fvalue, gdsq;

	fvalue = calc_f();
	gdsq = calc_g(gd);

	status = 0;
	pmf_orig = pmf;
	do {
		noNeg = true;
		for (int i = 1; i <= span; ++i) {
			pmf[i] = pmf_orig[i] + t * gd[i];
			noNeg = noNeg && (pmf[i] >= 0.0);
		}
		if (noNeg) {
			updateAux();
			new_fvalue = calc_f();
		}

		success = noNeg && new_fvalue >= fvalue + alpha * t * gdsq;
		if (status == 0) status = (success ? 2 : 1);
		
		if ((status == 1 && success) || (status == 2 && !success)) break;
		if (status == 1) t *= beta;
		else {
			t /= beta;
			if (t > max_step_size) break;
		}
		
	} while (true);
	
	if (!success) {
		t *= beta;
		for (int i = 1; i <= span; ++i) pmf[i] = pmf_orig[i] + t * gd[i];
	}
}

// Project pmf such that \sum pmf = 1
void FragLenDist::projection() {
	double sum = 0.0, value, minval = pmf[1];
	for (int i = 1; i <= span; ++i) {
		sum += pmf[i];
		if (minval > pmf[i]) minval = pmf[i];
	}
	
	if (fabs(sum - 1.0) < 1e-8) {
		for (int i = 1; i <= span; ++i) pmf[i] /= sum;
		return;
	}

	if (sum < 1.0) {
		value = (1.0 - sum) / span;
		for (int i = 1; i <= span; ++i) pmf[i] += value;
		return;
	}

	// sum > 1.0
	sum -= 1.0;
	value = sum / span;

	if (minval >= value) {
		for (int i = 1; i <= span; ++i) pmf[i] -= value;
		return;
	}

	for (int i = 1; i <= span; ++i) ids[i - 1] = i;
	sort(ids.begin(), ids.end(), *this);

	int i = 0, nleft = span;
	while (i < span && value > pmf[ids[i]]) {
		sum += pmf[ids[i]];
		pmf[ids[i]] = 0.0;
		++i; --nleft;
		assert(nleft > 0);
		value = sum / nleft;
	}
	for (; i < span; ++i) pmf[ids[i]] -= value;
}
