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

#include "Markov.hpp"
#include "Profile.hpp"
#include "QualDist.hpp"
#include "QProfile.hpp"
#include "SequencingErrorModel.hpp"

SequencingErrorModel::SequencingErrorModel(bool hasQual, bool to_log_space, int maxL) : hasQual(hasQual) {
	markov = NULL; profile = NULL; qprofile = NULL; qd = NULL;

	markov = new Markov(to_log_space);
	if (hasQual) qprofile = new QProfile(to_log_space);
	else profile = new Profile(to_log_space, maxL);
}

SequencingErrorModel::~SequencingErrorModel() {
	if (markov != NULL) delete markov;
	if (profile != NULL) delete profile;
	if (qprofile != NULL) delete qprofile;
}

void SequencingErrorModel::collect(const SequencingErrorModel* o) {
	markov->collect(o->markov);
	if (hasQual) qprofile->collect(o->qprofile);
	else profile->collect(o->profile);
}

void SequencingErrorModel::finish(int length) {
	markov->finish();
	if (hasQual) qprofile->finish();
	else profile->finish(length);
}

void SequencingErrorModel::clear() {
	markov->clear();
	if (hasQual) qprofile->clear();
	else profile->clear();
}

void SequencingErrorModel::read(std::ifstream& fin) {
	std::string line;

	// When RNASeqModel detects markov, it calls read of SequencingErrorModel
	markov->read(fin);
	
	while (getline(fin, line) && (line.length() == 0 || line[0] != '#'));

	if (hasQual) {
		assert(line.length() > 5 && line.substr(1, 4) == "qpro");
		qprofile->read(fin);		
	}
	else {
		assert(line.length() > 4 && line.substr(1, 3) == "pro");
		profile->read(fin);
	}
}

void SequencingErrorModel::write(std::ofstream& fout, bool isProb) {
	markov->write(fout, isProb);
	if (hasQual) qprofile->write(fout, isProb);
	else profile->write(fout, isProb);
}

void SequencingErrorModel::prepare_for_simulation(QualDist* qd) {
	this->qd = qd;
	markov->prepare_for_simulation();
	if (hasQual) qprofile->prepare_for_simulation();
	else profile->prepare_for_simulation();
}
