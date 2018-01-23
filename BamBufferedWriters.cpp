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

#include <cstdio>
#include <sstream>
#include <stdint.h>
#include <pthread.h>

#include "htslib/sam.h"
#include "utils.h"
#include "my_assert.h"
#include "BamWriter.hpp"
#include "BamBufferedWriters.hpp"

BamBufferedWriters::BamBufferedWriters(const char* imdName, int num_threads, int buffer_size) {
	char outF[STRLEN];

	// Allocate bam1_t* buffer
	this->buffer_size = buffer_size;
	buffer = new bam1_t*[buffer_size];
	for (int i = 0; i < buffer_size; ++i)
		buffer[i] = bam_init1();
	pos = (buffer_size > 1 ? 0 : -1); // if buffer_size == 1, pos = -1 to represent that we have no need to flush bam records
	
	// Allocate all BAM file writers  
	s = num_threads + 2;
	writers = new WriterType[s];

	// N0.bam
	sprintf(outF, "%s_N0.bam", imdName);
	writers[0].writer = new BamWriter(outF, NULL);
	writers[0].buffer = buffer;
	
	// thread_id.bam
	for (int i = 1; i <= num_threads; ++i) {
		sprintf(outF, "%s_%d.bam", imdName, i - 1);
		writers[i].writer = new BamWriter(outF, NULL);
		writers[i].buffer = buffer;
	}

	// N2.bam
	sprintf(outF, "%s_N2.bam", imdName);
	writers[s - 1].writer = new BamWriter(outF, NULL);
	writers[s - 1].buffer = buffer;

	// Set up pthreads
	threads.assign(s, pthread_t());
	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
}

BamBufferedWriters::~BamBufferedWriters() {
	if (buffer_size > 1) flush();
	else if (pos >= 0) writers[pos].write(buffer[0]);

	for (int i = 0; i < s; ++i) 
		delete writers[i].writer;	
	delete[] writers;
	for (int i = 0; i < buffer_size; ++i)
		bam_destroy1(buffer[i]);
	delete[] buffer;
	pthread_attr_destroy(&attr);  
}

void* write_to_bam_file(void* arg) {
	WriterType *writer = (WriterType*)arg;

	for (int i = 0; i < (int)writer->locations.size(); ++i)
		writer->write(writer->buffer[writer->locations[i]]);

	return NULL;
}

void BamBufferedWriters::flush() {
	if (pos == 0) return;

	// write BAM files
	for (int i = 0; i < s; ++i)
		if (writers[i].locations.size() > 0) {
			rc = pthread_create(&threads[i], &attr, write_to_bam_file, (void*)(&writers[i]));
			pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for writing BAM file!");
		}

	for (int i = 0; i < s; ++i)
		if (writers[i].locations.size() > 0) {
			rc = pthread_join(threads[i], NULL);
			pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for writing BAM file!");
		}

	pos = 0;
	for (int i = 0; i < s; ++i)
		writers[i].locations.clear();  
}
