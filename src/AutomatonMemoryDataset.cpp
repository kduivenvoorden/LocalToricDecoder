/*
 * AutomatonMemoryDataset.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "AutomatonMemoryDataset.h"

AutomatonMemoryDataset::AutomatonMemoryDataset(Point address, int Q, int U, int b)
	: countSignal(8,0), newCountSignal(8,0), flipSignal(8), newFlipSignal(9), count(9, std::vector<int>(b+1,0)),
	  address(address), Q(Q), U(U), b(b), isFlipedSignal(9,0),newIsFlipedSignal(8,0) {
	this->age = 0;
	for(int i = 0; i < 8; i++){
		this->flipSignal[i] = new Bit();
		this->newFlipSignal[i] = new Bit();
	}
	this->newFlipSignal[8] = new Bit();
}

AutomatonMemoryDataset::~AutomatonMemoryDataset() {
	for(int i = 0; i < 8;i++){
		delete flipSignal[i];
		delete newFlipSignal[i];
	}
}

