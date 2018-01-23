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

#ifndef BAMBUFFEREDWRITERS_H_
#define BAMBUFFEREDWRITERS_H_

#include <vector>
#include <pthread.h>

#include "htslib/sam.h"
#include "BamWriter.hpp"

struct WriterType {
	BamWriter *writer;
	bam1_t** buffer; // a pointer for BamBufferedWriter's buffer
	std::vector<int> locations;

	WriterType() {
		writer = NULL;
		buffer = NULL;
		locations.clear();
	}

	~WriterType() {
		if (writer != NULL) delete writer;
	}  

	void write(bam1_t* b) { writer->write(b); }
};

class BamBufferedWriters {
public:
	/*
		@param   imdName   intermediate file name
		@param   num_threads   number of threads, num_threads + 2 bam files will be generated
		@param   buffer_size   number of bam1_t alignments in the buffer
	 */
	BamBufferedWriters(const char* imdName, int num_threads, int buffer_size = 1000000);

	~BamBufferedWriters();

	bam1_t* get_bam1_t(int id) {
		if (buffer_size == 1) {
			if (pos >= 0) writers[pos].write(buffer[0]);
			pos = id;
			return buffer[0];
		}
		
		if (pos >= buffer_size) flush();
		writers[id].locations.push_back(pos);
		return buffer[pos++];
	}
	
private:
	int s, buffer_size, pos; // pos, next available buffer position, if buffer_size == 1, this is the last writer's id
	WriterType* writers;
	bam1_t** buffer;
	
	int rc;
	pthread_attr_t attr;
	std::vector<pthread_t> threads;

	void flush();
};

#endif
