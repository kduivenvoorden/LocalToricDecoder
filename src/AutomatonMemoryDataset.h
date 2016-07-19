/*
 * AutomatonMemoryDataset.h
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#ifndef AUTOMATONMEMORYDATASET_H_
#define AUTOMATONMEMORYDATASET_H_

#include "Bit.h"
#include <vector>

struct Point{
	Point(int x, int y) : x(x), y(y) {};
	int x;
	int y;
};

class AutomatonMemoryDataset {
public:
	int age;
	// we changed from arrays to vectors
	/*
	int* countSignal;
	int* newCountSignal;
	Bit** flipSignal;
	Bit** newFlipSignal;
	int** count;
	*/
	std::vector<int> countSignal;
	std::vector<int> newCountSignal;
	std::vector<Bit*> flipSignal;
	std::vector<Bit*> newFlipSignal;
	std::vector<std::vector<int> > count;
	Point address;
	int Q;
	int U;
	int b;

	AutomatonMemoryDataset(Point, int, int, int);
	virtual ~AutomatonMemoryDataset();
};

#endif /* AUTOMATONMEMORYDATASET_H_ */
