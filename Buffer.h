#ifndef BUFFER_H_
#define BUFFER_H_

#include<cstdio>
#include<fstream>

typedef unsigned long long bufsize_type;
const int FLOATSIZE = sizeof(float);

class Buffer {
public:
	Buffer(bufsize_type size, int sp, int nSamples, int cvlen, const char* tmpF) {
		cpos = 0;
		this->size = size;
		buffer = new float[size];
		ftmpOut.open(tmpF, std::ios::binary);

		fr = to = sp;
		this->nSamples = nSamples;
		this->cvlen = cvlen;

	}

	~Buffer() {
		if (fr < to) flushToTempFile();

		delete[] buffer;
		ftmpOut.close();
	}

	void write(float *vec) {
		if (size - cpos < bufsize_type(cvlen)) flushToTempFile();
		memcpy(buffer + cpos, vec, FLOATSIZE * cvlen);
		cpos += cvlen;
		++to;
	}

private:
	bufsize_type size, cpos; // cpos : current position

	float *buffer;
	std::ofstream ftmpOut;

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
