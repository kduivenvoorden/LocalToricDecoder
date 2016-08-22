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
	bool printMovie;
	std::vector<AutomatonMemoryDataset*> memory;
	Automaton** neighbors;
	Qubit** qubits;

	int localUpdateRules(Point, int[], Bit*[], int);
	int ipow(int, int);
	void printRule(int, int );
public:
	Automaton(Point, int, int, int, int, int, double, int, int,bool);
	virtual ~Automaton();

	void setNeighbors(Automaton**);
	void setQubits(Qubit**);
	void setSyndrome(int);

	void localDependencies();
	void executeFlips();
	void caUpdates();

};

#endif /* AUTOMATON_H_ */
