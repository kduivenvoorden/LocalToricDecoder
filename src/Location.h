/*
 * Location.h
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#ifndef LOCATION_H_
#define LOCATION_H_

class Location {
public:
	static const int N = 0;
	static const int W = 1;
	static const int E = 2;
	static const int S = 3;
	static const int NW = 4;
	static const int NE = 5;
	static const int SW = 6;
	static const int SE = 7;
	static const int C = 8;
	static int getOppositeDirection(int);

	static int oppositeDirections[9];
	Location();
	virtual ~Location();
};


#endif /* LOCATION_H_ */
