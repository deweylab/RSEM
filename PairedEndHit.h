#ifndef PAIREDENDHIT_H_
#define PAIREDENDHIT_H_

#include<iostream>

#include "SingleHit.h"

class PairedEndHit : public SingleHit {
public:
	PairedEndHit() : SingleHit() {
		insertL = 0;
	}

	PairedEndHit(int sid, int pos, int insertL, double conprb = 0.0) : SingleHit(sid, pos, conprb) {
		this->insertL = insertL;
	}

	int getInsertL() const { return insertL; }

	bool read(std::istream&);
	void write(std::ostream&);

private:
	int insertL; // insert length
};

bool PairedEndHit::read(std::istream& in) {
	conprb = 0.0;
    return (in>>sid>>pos>>insertL);
}

void PairedEndHit::write(std::ostream& out) {
	out<<" "<<sid<<" "<<pos<<" "<<insertL;
}

#endif /* PAIREDENDHIT_H_ */
