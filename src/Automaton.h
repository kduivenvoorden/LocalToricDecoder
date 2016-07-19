/*
 * Automaton.h
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#ifndef AUTOMATON_H_
#define AUTOMATON_H_

#include "AutomatonMemoryDataset.h"
#include "Qubit.h"
#include <vector>

class Automaton {
private:
	int hierarchyDepth;
	int syndromes[9];
	double fc;
	double fn;
	std::vector<AutomatonMemoryDataset*> memory;
	Automaton** neighbors;
	Qubit** qubits;

	void localUpdateRules(Point, int[], Bit*[], int);
	int ipow(int, int);
public:
	Automaton(Point, int, int, int, int, int, double, int, int);
	virtual ~Automaton();

	void setNeighbors(Automaton**);
	void setQubits(Qubit**);
	void setSyndrome(int);

	void localUpdates();
	void executeFlips();
	void signalUpdates();

};

#endif /* AUTOMATON_H_ */
