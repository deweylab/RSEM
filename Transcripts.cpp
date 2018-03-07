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

#include <string>
#include <fstream>

#include "utils.h"
#include "my_assert.h"

#include "Transcript.hpp"
#include "Transcripts.hpp"

void Transcripts::readFrom(const char* inpF) {
	std::string line;
	std::ifstream fin(inpF);

	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	fin>> M>> type;
	std::getline(fin, line);
	transcripts.assign(M + 1, Transcript());
	for (int i = 1; i <= M; ++i) {
		transcripts[i].read(fin);
	}
	fin.close();
}

void Transcripts::writeTo(const char* outF) {
	std::ofstream fout(outF);
	fout<< M<< " "<< type<< std::endl;
	for (int i = 1; i <= M; ++i) {
		transcripts[i].write(fout);
	}
	fout.close();
}

void Transcripts::writeTransListTo(const char* outF) {
	std::ofstream fout(outF);
	for (int i = 1; i <= M; ++i)
		fout<< transcripts[i].getTranscriptID()<< '\t'<< transcripts[i].getLength()<< std::endl;
	fout.close();
}

void Transcripts::updateCLens() {
	for (int i = 1; i <= M; ++i)
		transcripts[i].updateCLen();
}
