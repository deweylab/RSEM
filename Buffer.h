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
	Buffer(int nMB, int nSamples, int cvlen, const char* tmpF) {
		cpos = 0;
		size = bufsize_type(nMB) * 1024 * 1024 / FLOATSIZE / cvlen;
		if (size > (bufsize_type)nSamples) size = nSamples;
		general_assert(size > 0, "Memory allocated for credibility intervals is not enough!");
		size *= cvlen;

		buffer = new float[size];
		ftmpOut.open(tmpF, std::ios::binary);
		pthread_mutex_init(&lock, NULL);

		fr = to = 0;
		this->nSamples = nSamples;
		this->cvlen = cvlen;
	}

	~Buffer() {
		if (fr < to) flushToTempFile();

		delete[] buffer;
		pthread_mutex_destroy(&lock);
		ftmpOut.close();
	}

	void write(int n, float **vecs) {
		pthread_assert(pthread_mutex_lock(&lock), "pthread_mutex_lock", "Error occurred while acquiring the lock!");
		for (int i = 0; i < n; i++) {
			if (size - cpos < bufsize_type(cvlen)) flushToTempFile();
			memcpy(buffer + cpos, vecs[i], FLOATSIZE * cvlen);
			cpos += cvlen;
			++to;
		}
		pthread_assert(pthread_mutex_unlock(&lock), "pthread_mutex_unlock", "Error occurred while releasing the lock!");
	}

private:
	bufsize_type size, cpos; // cpos : current position

	float *buffer;
	std::ofstream ftmpOut;
	pthread_mutex_t lock;

	int fr, to; // each flush, sample fr .. to - 1
	int nSamples, cvlen;

	void flushToTempFile() {
		std::streampos gap1 = std::streampos(fr) * FLOATSIZE;
		std::streampos gap2 = std::streampos(nSamples - to) * FLOATSIZE;
		float *p = NULL;

		ftmpOut.seekp(0, std::ios::beg);
		for (int i = 0; i < cvlen; i++) {
			p = buffer + i;
			ftmpOut.seekp(gap1, std::ios::cur);
			for (int j = fr; j < to; j++) {
				ftmpOut.write((char*)p, FLOATSIZE);
				p += cvlen;
			}
			ftmpOut.seekp(gap2, std::ios::cur);
		}

		cpos = 0;
		fr = to;
	}
};

#endif /* BUFFER_H_ */
