/*
 * Simulator.h
 *
 *  Created on: Aug 20, 2015
 *      Author: michels
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include "ToricCode.h"
#include "Automaton.h"
#include <vector>
#include <list>

class Simulator {
private:
	int L;
	int Q;
	int k;
	int U;
	int b;
	int tmax;
	int N;
	double p;
	double q;
	int fn;
	int fc;
	ToricCode* code;
	Automaton*** automata;

	std::vector<double> linspace(double, double, int);
	std::vector<int> linspace(int, int, int);

	void init();
	void reset();
	void simulationStep();
	/* simulation critical methods are also implemented to work without the member variables
	 * code and automata. this is neccessary for using openmp because it is not possible
	 * to give every thread a private copy of a member variable.
	 */
	void simulationStep(int**);
	void reset(ToricCode*&, Automaton***&, int);
	void init(ToricCode*&, Automaton***&, int, int, int, int, int, double, double, int, int);
	void simulationStep(ToricCode*&, Automaton***&, int);

public:
	/* functions implementing tests */
	void level0ErrorCorrectionTest();
	void level1ErrorCorrectionTest();
	void level0SyndromeCorrectionTest();
	/* funtion to generate a threshold curve for the harrington decoder */
	std::list<double> generateThresholdCurve(int, int, int, int, int, int, double, double, int, int);
	/* funtion to generate a threshold curve for the minimum weight perfect matching decoder */
	void generateMWMThresholdCurve(int, int, int, double, double);
	/* funtion to generate a decay time curve for the harrington decoder */
	std::list<int> generateDecayCurve(int, int, int, int, int, int, double, int, int);
	/* funtion to generate a decay time curve for the minimum weight perfect matching decoder */
	std::list<int> generateMWMDecayCurve(int, double);

	int getTmax();
	int getN();

	void setTmax(int);
	void setN(int);

	Simulator();
	virtual ~Simulator();
};

#endif /* SIMULATOR_H_ */
