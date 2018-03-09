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

#include <new>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "MateLenDist.hpp"
#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"
#include "QualDist.hpp"
#include "NoiseProfile.hpp"
#include "NoiseQProfile.hpp"
#include "IlluminaSequenceModel.hpp"

IlluminaSequenceModel::IlluminaSequenceModel(model_mode_type mode, bool has_qual, int min_len, int max_len) : mode(mode), has_qual(has_qual) {
	mld = NULL; markov = NULL; pro = NULL; qpro = NULL; qd = NULL; npro = NULL; nqpro = NULL;

	if (mode == FIRST_PASS || mode == MASTER || mode == SIMULATION) mld = new MateLenDist(mode, min_len, max_len);
	markov = new Markov(mode);
	if (has_qual) {
		if (mode == FIRST_PASS || mode == MASTER || mode == SIMULATION) qd = new QualDist(mode);
		qpro = new QProfile(mode);
		nqpro = new NoiseQProfile(mode);
	}
	else {
		pro = new Profile(mode, max_len);
		npro = new NoiseProfile(mode, max_len);
	}
}

IlluminaSequenceModel::~IlluminaSequenceModel() {
	if (mld != NULL) delete mld;
	if (markov != NULL) delete markov;
	if (qd != NULL) delete qd;
	if (qpro != NULL) delete qpro;
	if (nqpro != NULL) delete nqpro;
	if (pro != NULL) delete pro;
	if (npro != NULL) delete npro;
}

void IlluminaSequenceModel::clear() {
	if (mode == FIRST_PASS) mld->clear();
	markov->clear();
	if (has_qual) {
		if (mode == FIRST_PASS) qd->clear();
		qpro->clear();
		nqpro->clear();
	}
	else {
		pro->clear();
		npro->clear();
	}
}

void IlluminaSequenceModel::collect(const IlluminaSequenceModel* o) {
	markov->collect(o->markov);
	if (has_qual) {
		qpro->collect(o->qpro);
		nqpro->collect(o->nqpro);
	}
	else {
		pro->collect(o->pro);
		npro->collect(o->npro);
	}
}

void IlluminaSequenceModel::finish() {
	if (mode == FIRST_PASS) {
		mld->findBoundaries();
		if (!has_qual) {
			pro->setMaxL(mld->getMaxL());
			npro->setMaxL(mld->getMaxL());
		}
		mld->finish();
	}
	markov->finish();
	if (has_qual) {
		if (mode == FIRST_PASS) qd->finish();
		qpro->finish();
		nqpro->finish();
	}
	else {
		pro->finish();
		npro->finish();
	}
}

void IlluminaSequenceModel::read(std::ifstream& fin, int choice) {
	std::string line;

	assert((fin>> line) && (line == "#IlluminaSequenceModel:"));
	getline(fin, line);

	if (choice < 2) {
		mld->read(fin, choice);
		markov->read(fin, choice);
		if (has_qual) {
			qd->read(fin, choice);
			qpro->read(fin, choice);
		}
		else pro->read(fin, choice);
	}

	if (has_qual) nqpro->read(fin, choice);
	else npro->read(fin, choice);
}

void IlluminaSequenceModel::write(std::ofstream& fout, int choice) {
	fout<< "#IlluminaSequenceModel:";
	if (choice < 2) fout<< "\tMateLenDist\tMarkov"<< (has_qual ? "\tQualDist\tQProfile" : "\tProfile");
	fout<< (has_qual ? "\tNoiseQProfile" : "\tNoiseProfile")<< std::endl<< std::endl;
	
	if (choice < 2) {
		mld->write(fout, choice);
		markov->write(fout, choice);
		if (has_qual) {
			qd->write(fout, choice);
			qpro->write(fout, choice);
		}
		else pro->write(fout, choice);
	}

	if (has_qual) nqpro->write(fout, choice);
	else npro->write(fout, choice);
}
