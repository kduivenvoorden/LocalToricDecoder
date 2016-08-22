/*
 * ToricCode.h
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#ifndef TORICCODE_H_
#define TORICCODE_H_

#include "Qubit.h"
#include "Location.h"
#include <random>

class ToricCode {
private:
	int L;
	int** vertices;
	Qubit*** edges;
	int** starSyndromes;
	int** faultySyndromes;
	double p;
	double q;
	bool printMovie;

	std::random_device generator;
	std::mt19937 mt;
	std::uniform_real_distribution<double> distribution;

	bool hasLogicalError(Qubit***, int);
	void flip(Qubit***, int, int, int);
	void dump(int**, Qubit***, int);
	Qubit* getQubit(Qubit***, int, int, int, int);
	
	
public:
	ToricCode(int, double, double, bool);
	virtual ~ToricCode();
	void reset();
	int index(int);
	int starCheck(int, int);
	void flip(int, int, int);
	int** getStarSyndromes();
	void introduceErrors();
	bool hasSyndrome();
	bool hasLogicalError();
	bool hasLogicalErrorAfterMWM();
	void minimumWeightMatching();
	Qubit** getQubits(int, int);
	void dump(void);
	int getDistance(int,int,int,int);
	void setQ(double);
	void viewQubitErrors();
	void viewSyndromeErrors();
	void viewStarSyndromes();
	void viewRealStarSyndromes();
};

#endif /* TORICCODE_H_ */
