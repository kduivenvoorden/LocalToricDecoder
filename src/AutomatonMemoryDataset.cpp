/*
 * AutomatonMemoryDataset.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "AutomatonMemoryDataset.h"

AutomatonMemoryDataset::AutomatonMemoryDataset(Point address, int Q, int U, int b)
	: countSignal(8,0), newCountSignal(8,0), flipSignal(8), newFlipSignal(8), count(9, std::vector<int>(9,0)),
	  address(address), Q(Q), U(U), b(b) {
	this->age = 0;
	for(int i = 0; i < 8; i++){
		this->flipSignal[i] = new Bit();
		this->newFlipSignal[i] = new Bit();
	}
	/*this->countSignal = new int[8];
	this->newCountSignal = new int[8];
	this->flipSignal = new Bit*[8];
	this->newFlipSignal = new Bit*[8];
	this->count = new int*[9];
	for(int i = 0; i < 9; i++){
		this->count[i] = new int[2];
	}
	*/
}

AutomatonMemoryDataset::~AutomatonMemoryDataset() {
	/*delete[] countSignal;
	delete[] newCountSignal;
	for(int i = 0; i< 8; i++)
		delete flipSignal[i];
	delete[] flipSignal;
	for(int i = 0; i< 8; i++)
		delete newFlipSignal[i];
	delete[] newFlipSignal;
	for(int i = 0; i < 9; i++)
		delete[] count[i];
	delete[] count;*/
	for(int i = 0; i < 8;i++){
		delete flipSignal[i];
		delete newFlipSignal[i];
	}
}

