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

#include <cassert>
#include <string>
#include <fstream>

#include "MateLenDist.hpp"
#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"
#include "QualDist.hpp"
#include "NoiseProfile.hpp"
#include "NoiseQProfile.hpp"
#include "IlluminaSequenceModel.hpp"

IlluminaSequenceModel::IlluminaSequenceModel(int mode, int status, int min_len, int max_len) : mode(mode), status(status), min_len(min_len), max_len(max_len) {
	mld = NULL; markov = NULL; pro = NULL; qpro = NULL; qd = NULL; npro = NULL; nqpro = NULL;

	if (mode != 1) mld = new MateLenDist(mode, min_len, max_len);
	markov = new Markov(mode);
	if (status & 1) {
		if (mode != 1) qd = new QualDist(mode);
		qpro = new QProfile(mode);
		nqpro = new NoiseQProfile(mode);
	}
	else {
		pro = new Profile(mode, maxL);
		npro = new NoiseProfile(mode, maxL);
	}
}

IlluminaSequenceModel::~IlluminaSequenceModel() {
	if (mode != 1) delete mld;
	delete markov;
	if (status & 1) {
		if (mode != 1) delete qd;
		delete qpro;
		delete nqpro;
	}
	else {
		delete pro;
		delete npro;
	}
}

void IlluminaSequenceModel::clear() {
	if (mode == 0 && (status & 2)) mld->clear();
	markov->clear();
	if (status & 1) {
		if (mode == 0 && (status & 2)) qd->clear();
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
	if (status & 1) {
		qpro->collect(o->qpro);
		nqpro->collect(o->nqpro);
	}
	else {
		pro->collect(o->pro);
		npro->collect(o->npro);
	}
}

void IlluminaSequenceModel::finish() {
	if (mode == 0 && (status & 2)) mld->finish();
	markov->finish();
	if (status & 1) {
		if (mode == 0 && (status & 2)) qd->finish();
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

	fin>> line;
	assert(line == "#IlluminaSequenceModel:");
	getline(fin, line);

	if (choice < 2) {
		mld->read(fin, choice);
		markov->read(fin, choice);
		if (status & 1) {
			qd->read(fin, choice);
			qpro->read(fin, choice);
		}
		else pro->read(fin, choice);
	}

	if (status & 1) nqpro->read(fin, choice);
	else npro->read(fin, choice);
}

void IlluminaSequenceModel::write(std::ofstream& fout, int choice) {
	fout<< "#IlluminaSequenceModel:";

	if (choice < 2) {
		mld->write(fout, choice); fout<< "\tMateLenDist";
		markov->write(fout, choice); fout<< "\tMarkov";
		if (status & 1) {
			qd->write(fout, choice); fout<< "\tQualDist";
			qpro->write(fout, choice); fout<< "\tQProfile";
		}
		else {
			pro->write(fout, choice); fout<< "\tProfile";
		}
	}

	if (status & 1) { nqpro->write(fout, choice); fout<< "\tNoiseQProfile"; }
	else { npro->write(fout, choice); fout<< "\tNoiseProfile"; }

	fout<< std::endl;
}
