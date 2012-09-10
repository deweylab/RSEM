#ifndef BUFFER_H_
#define BUFFER_H_

#include<cstdio>
#include<fstream>
#include<pthread.h>

#include "my_assert.h"

typedef unsigned long long bufsize_type;
const int FLOATSIZE = sizeof(float);

class Buffer {
public:
	// in_mem_arr must be allocated memory before the Buffer is constructed
	Buffer(int nMB, int nSamples, int vlen, float* in_mem_arr, const char* tmpF) {
		cpos = 0;
		size = bufsize_type(nMB) * 1024 * 1024 / FLOATSIZE / vlen;
		if (size > (bufsize_type)nSamples) size = nSamples;
		general_assert(size > 0, "Memory allocated for credibility intervals is not enough!");
		size *= vlen;

		buffer = new float[size];
		ftmpOut.open(tmpF, std::ios::binary);
		pthread_mutex_init(&lock, NULL);

		fr = to = 0;
		this->nSamples = nSamples;
		this->vlen = vlen;
		this->in_mem_arr = in_mem_arr;
	}

	~Buffer() {
		if (fr < to) flushToTempFile();

		delete[] buffer;
		pthread_mutex_destroy(&lock);
		ftmpOut.close();
	}

	void write(float value, float *vec) {
		pthread_assert(pthread_mutex_lock(&lock), "pthread_mutex_lock", "Error occurred while acquiring the lock!");
		if (size - cpos < bufsize_type(vlen)) flushToTempFile();
		in_mem_arr[to] = value;
		memcpy(buffer + cpos, vec, FLOATSIZE * vlen);
		cpos += vlen;
		++to;
		pthread_assert(pthread_mutex_unlock(&lock), "pthread_mutex_unlock", "Error occurred while releasing the lock!");
	}

private:
	bufsize_type size, cpos; // cpos : current position

	float *buffer;
	float *in_mem_arr;
	std::ofstream ftmpOut;
	pthread_mutex_t lock;

	int fr, to; // each flush, sample fr .. to - 1
	int nSamples, vlen; // vlen : vector length

	void flushToTempFile() {
		std::streampos gap1 = std::streampos(fr) * FLOATSIZE;
		std::streampos gap2 = std::streampos(nSamples - to) * FLOATSIZE;
		float *p = NULL;

		ftmpOut.seekp(0, std::ios::beg);
		for (int i = 0; i < vlen; i++) {
			p = buffer + i;
			ftmpOut.seekp(gap1, std::ios::cur);
			for (int j = fr; j < to; j++) {
				ftmpOut.write((char*)p, FLOATSIZE);
				p += vlen;
			}
			ftmpOut.seekp(gap2, std::ios::cur);
		}

		cpos = 0;
		fr = to;
	}
};

#endif /* BUFFER_H_ */
