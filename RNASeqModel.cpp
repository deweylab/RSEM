/* Copyright (c) 2016
   Bo Li (University of California, Berkeley)
   bli25@berkeley.edu

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

#include<fstream>

#include "my_assert.h"
#include "RNASeqModel.hpp"

RNASeqModel::RNASeqModel(int model_type) : model_type(model_type) {
	mld1 = mld2 = NULL;
	npro = NULL; qd = NULL; nqpro = NULL;
	frag_min = 99999999; frag_max = 0;
	
	mld1 = new MateLenDist();
	if (model_type >= 2) mld2 = new MateLenDist();
	if (model_type & 1) {
		qd = new QualDist(true); // enable sufficient statistics
		nqpro = new NoiseQProfile();
	}
	else npro = new NoiseProfile();
}

RNASeqModel::~RNASeqModel() {
	if (mld1 != NULL) delete mld1;
	if (mld2 != NULL) delete mld2;
	if (npro != NULL) delete npro;
	if (qd != NULL) delete qd;
	if (nqpro != NULL) delete nqpro;
}

// ssF: model_type \n mld1 \n mld2 \n qd \n nqpro \n pro
// paramF : frag_min frag_max \n maxL
void RNASeqModel::write(const char* ssF, const char* paramF) {
	// ssF
	std::ofstream fout(ssF);
	assert(fout.is_open());

	fout<< model_type<< std::endl<< std::endl;
	mld1->write(fout, false);
	if (model_type >= 2) mld2->write(fout, false);
	if (model_type & 1) {
		qd->finish();
		qd->write(fout);
		nqpro->write(fout);
	}
	else npro->write(fout);

	fout.close();

	// paramF
	fout.open(paramF);
	assert(fout.is_open());
	fout<< frag_min<< '\t'<< frag_max<< std::endl;
	int maxL = mld1->getMaxL();
	if (model_type >= 2) {
		general_assert(mld2->getMinL() <= mld2->getMaxL(), "Model type is " + itos(model_type) + ", but no paired-end reads are detected!");
		if (maxL < mld2->getMaxL()) maxL = mld2->getMaxL();
	}
	fout<< maxL<< std::endl;
	fout.close();
}

