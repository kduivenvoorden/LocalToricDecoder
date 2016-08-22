/*
 * Simulator.cpp
 *
 *  Created on: Aug 20, 2015
 *      Author: michels
 */

#include "Simulator.h"
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <vector>

Simulator::Simulator() {}
Simulator::~Simulator() {}

int Simulator::getN(){
	return this->N;
}

int Simulator::getTmax(){
	return this->tmax;
}

void Simulator::setTmax(int tmax){
	this->tmax = tmax;
}

void Simulator::setN(int N){
	this->N = N;
}

void Simulator::setPrintMovie(bool printMovie){
	this->printMovie = printMovie;
}

std::list<int> Simulator::generateDecayCurve(int L, int Q, int k, int U, int b, int tau, double p, int fn, int fc, bool printMovie){
	this->L = L;
	this->Q = Q;
	this->U = U;
	this->b = b;
	this->k = k;
	this->fn = fn;
	this->fc = fc;
	this->printMovie = printMovie;

	this->p = p;
	double q = p;

	int tmax = this->tmax;
	int N = this->N;
	bool doNotRoughTest = false;

	std::list<int> decayTimes;
	//long sum = 0;
	//int globalIndex = 1;

	ToricCode* code;
	Automaton*** automata;
	//#pragma omp parallel for private(code,automata)
	for(int trial = 0; trial < N; trial++){
		//if(trial % 100 ==0){
		//std::cout << "Running Trial " << trial << " out of " << N << std::endl;
		//}
		this->init(code,automata,L,Q,U,b,k,p,q,fn,fc,printMovie);
		if (printMovie)
			std::cout << "-1" << std::endl;
		for(int time = 0; time < tmax; time++){
			if(time % tau == 0){
				if (printMovie)
					code-> viewQubitErrors();
				code-> introduceErrors();
				if (printMovie) {
					code-> viewQubitErrors();
					code-> viewSyndromeErrors();
				}
			} else {
				if (printMovie) {
					code-> viewQubitErrors();				
					code-> viewQubitErrors();
					code-> viewSyndromeErrors();
				}
			}
			code->getStarSyndromes();
			if (printMovie)
				code-> viewStarSyndromes();
			
			if(!code->hasSyndrome()){
				time = time + (tau - time % tau-1);
				if (printMovie)
					std::cout << "-1" << std::endl;
				continue;
			}
			this->simulationStep(code,automata,L);
			if (printMovie) 
				std::cout << "-1" << std::endl;
			
			if (doNotRoughTest || code->hasLogicalError()) {
				code->getStarSyndromes();
				if(code->hasLogicalErrorAfterMWM()){
					//#pragma omp critical
					{
						decayTimes.push_back(time/tau);
						if(!printMovie)
							std::cout << time/tau << std::endl;
						//sum += time;
						//std::cout << "  Average so far: " << sum/globalIndex << std::endl;
						//globalIndex++;
						//break;
						time = tmax;
					}
				}
			}
		}
		if (printMovie)
			std::cout << "-1" << std::endl;
		this->reset(code, automata, L);
	}
	return decayTimes;
}

void Simulator::init() {
	this->init(this->code, this->automata, this->L, this->Q, this->U, this->b, this->k,
			this->p, this->q, this->fn, this->fc, this->printMovie);
}

void Simulator::init(ToricCode* &code, Automaton*** &automata, int L, int Q, int U, int b, int k,
		double p, double q, int fn, int fc,bool printMovie){
	code = new ToricCode(L, p, q, printMovie);
	automata = new Automaton**[L];
	for (int i = 0; i < L; i++) {
		automata[i] = new Automaton*[L];
		for (int j = 0; j < L; j++) {
			automata[i][j] = new Automaton(Point(j, L - i - 1), L, Q, U, b, k,p,fn,fc,printMovie);
		}
	}

	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			Automaton** neighbors = new Automaton*[8];
			neighbors[0] = automata[code->index(i - 1)][code->index(j)];
			neighbors[1] = automata[code->index(i)][code->index(j - 1)];
			neighbors[2] = automata[code->index(i)][code->index(j + 1)];
			neighbors[3] = automata[code->index(i + 1)][code->index(j)];
			neighbors[4] = automata[code->index(i - 1)][code->index(j - 1)];
			neighbors[5] = automata[code->index(i - 1)][code->index(j + 1)];
			neighbors[6] = automata[code->index(i + 1)][code->index(j - 1)];
			neighbors[7] = automata[code->index(i + 1)][code->index(j + 1)];

			automata[i][j]->setNeighbors(neighbors);
			automata[i][j]->setQubits(code->getQubits(i, j));
		}
	}
}

void Simulator::reset() {
	this->reset(this->code, this->automata, this->L);
}

void Simulator::reset(ToricCode* &code, Automaton*** &automata, int L) {
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			delete automata[i][j];
		}
		delete[] automata[i];
	}
	delete[] automata;
	delete code;
}

void Simulator::simulationStep(ToricCode*& code, Automaton***& automata, int L) {
	int** syndromes = code->getStarSyndromes();
	
		
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->setSyndrome(syndromes[i][j]);
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->localDependencies();
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->executeFlips();
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->caUpdates();
		}
	}

	
}

void Simulator::simulationStep() {
	int** syndromes = code->getStarSyndromes();
	this->simulationStep(syndromes);
}

void Simulator::simulationStep(int **syndromes) {
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->setSyndrome(syndromes[i][j]);
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->localDependencies();
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->executeFlips();
		}
	}
	
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->caUpdates();
		}
	}

}

std::vector<double> Simulator::linspace(double a, double b, int n){
	std::vector<double> array(n);
	double step = (b-a) / (n-1);

	for(int i = 0; i < n; i++){
	    array[i] = a+i*step;
	}
	return array;
}

std::vector<int> Simulator::linspace(int a, int b, int n){
	std::vector<int> array(n);
	double step = (b-a) / (n-1);

	for(int i = 0; i < n; i++){
	    array[i] = (a+i*step);
	}
	return array;
}
