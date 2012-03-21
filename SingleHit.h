#ifndef SINGLEHIT_H_
#define SINGLEHIT_H_

#include<cstdlib>
#include<iostream>

//char dir : 0 +, 1 - , encoding as 1 + , -1 -
class SingleHit {
public:
	SingleHit() {
		sid = 0; pos = -1; conprb = 0.0; // for noise gene
	}

	//sid encodes dir here
	SingleHit(int sid, int pos, double conprb = 0.0) {
		this->sid = sid;
		this->pos = pos;
		this->conprb = conprb;
	}

	bool isNoise() const { return sid == 0; }

	//makes no sense for noise gene
	int getDir() const { return sid < 0; }

	int getSid() const { return abs(sid); }

	int getPos() const { return pos; }

	double getConPrb() const { return conprb; }

	void setConPrb(double conprb) {
	    this->conprb = conprb;
	}

	bool read(std::istream&);
	void write(std::ostream&);

protected:
	int sid, pos; // sid encodes dir
	double conprb; // conditional probability
};

bool SingleHit::read(std::istream& in) {
	conprb = 0.0;
	return (in>>sid>>pos);
}

void SingleHit::write(std::ostream& out) {
	out<<" "<<sid<<" "<<pos;
}

#endif /* SINGLEHIT_H_ */



