/*
 * Location.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "Location.h"



Location::Location() {
		// TODO Auto-generated constructor stub
}

Location::~Location() {
	// TODO Auto-generated destructor stub
}

int Location::getOppositeDirection(int direction){
	/*switch(direction){
	case Location::N:
		return Location::S;
	case Location::S:
		return Location::N;
	case Location::W:
		return Location::E;
	case Location::E:
		return Location::W;
	case Location::NW:
		return Location::SE;
	case Location::SE:
		return Location::NW;
	case Location::NE:
		return Location::SW;
	case Location::SW:
		return Location::NE;
	default:
		return Location::C;
	}
	*/

	return Location::oppositeDirections[direction];
}

int Location::oppositeDirections[] = {Location::S, Location::E, Location::W, Location::N, Location::SE, Location::SW, Location::NE, Location::NW, Location::C};
