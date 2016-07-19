/*
 * Bit.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "Bit.h"

void Bit::flip(){
	this->bit = (this->bit + 1) % 2;
}

Bit::Bit() {
	this->bit = 0;
}

Bit::~Bit() {
}

