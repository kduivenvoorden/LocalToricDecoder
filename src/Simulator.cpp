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

void Simulator::level0ErrorCorrectionTest(){
	this->L = 4;
	this->Q = 4;
	this->k = 1;
	this->U = 4;
	this->b = 2;
	this->p = 0.00;
	this->q = this->p;

	this->init();

	auto step = [this](){
		std::cout << "##########################################" << std::endl;
		this->code->getStarSyndromes();
		this->code->dump();
		std::cout << std::endl;
		for(int i = 0; i < 3; i++){
			this->simulationStep();
			this->code->getStarSyndromes();
			this->code->dump();
			std::cout << std::endl;
		}
	};

	this->code->flip(1,1,Location::S);
	step();

	this->code->flip(1,1,Location::E);
	step();

	this->code->flip(1,1,Location::E);
	this->code->flip(1,1,Location::S);
	step();

	this->code->flip(1,1,Location::S);
	this->code->flip(1,1,Location::E);
	step();

	this->code->flip(1,1,Location::S);
	this->code->flip(1,2,Location::S);
	step();

	this->code->flip(1,1,Location::E);
	this->code->flip(2,1,Location::E);
	step();

	this->reset();
}

void Simulator::level1ErrorCorrectionTest(){
	this->L = 16;
	this->Q = 4;
	this->k = 2;
	this->U = 4;
	this->b = 2;
	this->p = 0.00;
	this->q = this->p;

	this->init();

	auto step = [this](){
		std::cout << "##########################################" << std::endl;
		this->code->getStarSyndromes();
		this->code->dump();
		std::cout << std::endl;
		for(int i = 0; i < 8*4; i++){
			this->simulationStep();
			std::cout << "After Day: " << i << std::endl;
			this->code->getStarSyndromes();
			this->code->dump();
			std::cout << std::endl;
		}
	};

	this->code->flip(5,6,Location::S);
	this->code->flip(6,6,Location::S);
	this->code->flip(7,6,Location::S);
	this->code->flip(8,6,Location::S);
	step();

	this->code->flip(5,6,Location::E);
	this->code->flip(5,7,Location::E);
	this->code->flip(5,8,Location::E);
	this->code->flip(5,9,Location::E);
	step();

	this->code->flip(5,6,Location::E);
	this->code->flip(5,7,Location::E);
	this->code->flip(5,8,Location::E);
	this->code->flip(5,9,Location::E);
	this->code->flip(5,6,Location::S);
	this->code->flip(6,6,Location::S);
	this->code->flip(7,6,Location::S);
	this->code->flip(8,6,Location::S);
	step();

	this->code->flip(9,6,Location::E);
	this->code->flip(9,7,Location::E);
	this->code->flip(9,8,Location::E);
	this->code->flip(9,9,Location::E);
	this->code->flip(5,6,Location::S);
	this->code->flip(6,6,Location::S);
	this->code->flip(7,6,Location::S);
	this->code->flip(8,6,Location::S);
	step();

	this->code->flip(9,6,Location::E);
	this->code->flip(9,7,Location::E);
	this->code->flip(9,8,Location::E);
	this->code->flip(9,9,Location::E);
	this->code->flip(5,6,Location::E);
	this->code->flip(5,7,Location::E);
	this->code->flip(5,8,Location::E);
	this->code->flip(5,9,Location::E);
	step();

	this->code->flip(5,6,Location::S);
	this->code->flip(6,6,Location::S);
	this->code->flip(7,6,Location::S);
	this->code->flip(8,6,Location::S);
	this->code->flip(5,10,Location::S);
	this->code->flip(6,10,Location::S);
	this->code->flip(7,10,Location::S);
	this->code->flip(8,10,Location::S);
	step();

	this->reset();

}

void Simulator::level0SyndromeCorrectionTest(){
	this->L = 4;
	this->Q = 4;
	this->k = 1;
	this->U = 4;
	this->b = 2;
	this->p = 0.00;
	this->q = this->p;

	this->init();

	int** syndromes;

	auto step = [this](int** syndromes){
		for(int i = 0; i < 3; i++){
			this->simulationStep(syndromes);
			syndromes = this->code->getStarSyndromes();
			this->code->dump();
			std::cout << std::endl;
		}
	};

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[2][1] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[1][2] = 1;
	syndromes[1][1] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[1][1] = 1;
	syndromes[2][2] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[2][1] = 1;
	syndromes[1][2] = 1;
	syndromes[2][2] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[1][1] = 1;
	syndromes[1][2] = 1;
	syndromes[2][2] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	std::cout << "##########################################" << std::endl;
	syndromes = this->code->getStarSyndromes();
	syndromes[1][1] = 1;
	syndromes[1][2] = 1;
	syndromes[2][1] = 1;
	syndromes[2][2] = 1;
	code->dump();
	std::cout << std::endl;
	step(syndromes);

	this->reset();
}

void Simulator::generateMWMThresholdCurve(int L, int Q, int resolution, double pmin, double pmax){
	this->L = L;
	this->Q = Q;

	/*int resolution = 20;
	double pmin = 0.01;
	double pmax = 0.04;*/
	std::vector<double> ps = this->linspace(pmin,pmax,resolution);

	for(int i = 0; i < resolution; i++){
		p = ps[i];
		this->q = 0;

		this->code = new ToricCode(L,p,this->q);
		int numberOfFailures = 0;
		for(int trial = 0; trial < this->N; trial++){
			this->code->introduceErrors();
			this->code->getStarSyndromes();
			this->code->minimumWeightMatching();
			this->code->getStarSyndromes();
			if(this->code->hasLogicalError()){
				numberOfFailures++;
			}
			this->code->reset();
		}
		delete this->code;
		std::cout << std::setw(12) << p << ":" << std::setw(12) << numberOfFailures << std::endl;
	}
}

std::list<double> Simulator::generateThresholdCurve(int L, int Q, int k, int U, int b, int resolution, double pmin, double pmax,
	int fn, int fc){
	this->L = L;
	this->Q = Q;
	this->k = k;
	this->U = U;
	this->b = b;
	this->fn = fn;
	this->fc = fc;

	/*int resolution = 20;
	double pmin = 0.01;
	double pmax = 0.04;*/
	std::vector<double> ps = this->linspace(pmin,pmax,resolution);

	std::list<double> plogicalList;
	for(int i = 0; i < resolution; i++){
		double p = ps[i];
		double q = 0;

		int numberOfFailures = 0;
		ToricCode* code;
		Automaton*** automata;
		#pragma omp parallel for reduction(+:numberOfFailures) private(code,automata)
		for(int trial = 0; trial < this->N; trial++){
			this->init(code,automata,L,Q,U,b,k,p,q,fn,fc);
			code->introduceErrors();
			for(int time = 0; time < this->tmax; time++){
				this->simulationStep(code,automata,L);
				code->getStarSyndromes();
				if(!code->hasSyndrome()){
					break;
				}
			}
			if(code->hasLogicalError()){
				numberOfFailures++;
			}
			this->reset(code,automata,L);
		}
		std::cout << p << "," << (double)numberOfFailures/((double)this->N)<< std::endl;
	}
	return plogicalList;
}

std::list<int> Simulator::generateMWMDecayCurve(int L, double p){
	this->L = L;
	this->q = 0;
	this->p = p;

	int tmax = this->tmax;
	int N = this->N;


	std::list<int> decayTimes;
	long sum = 0;
	int globalIndex = 1;
	for(int trial = 0; trial < N; trial++){
		for(int time = 0; time < tmax; time++){
			this->code = new ToricCode(L, this->p, this->q);
			this->code->introduceErrors();
			this->code->getStarSyndromes();
			this->code->minimumWeightMatching();
			this->code->getStarSyndromes();
			if(this->code->hasSyndrome()){
				std::cout << time << std::endl;
				this->code->dump();
				return decayTimes;
			}
			if(this->code->hasLogicalError()){
				decayTimes.push_back(time);
				std::cout << "Decay time: " << time << std::endl;
				sum += time;
				std::cout << "Average so far: " << sum/globalIndex << std::endl;
				globalIndex++;
				break;
			}
			delete this->code;
		}
	}
	return decayTimes;
}

std::list<int> Simulator::generateDecayCurve(int L, int Q, int k, int U, int b, int tau, double p, int fn, int fc){
	this->L = L;
	this->Q = Q;
	this->U = U;
	this->b = b;
	this->k = k;
	this->fn = fn;
	this->fc = fc;

	this->p = p;
	double q = p;

	int tmax = this->tmax;
	int N = this->N;

	std::list<int> decayTimes;
	//long sum = 0;
	//int globalIndex = 1;

	ToricCode* code;
	Automaton*** automata;
	#pragma omp parallel for private(code,automata)
	for(int trial = 0; trial < N; trial++){
		//if(trial % 100 ==0){
		//std::cout << "Running Trial " << trial << " out of " << N << std::endl;
		//}
		this->init(code,automata,L,Q,U,b,k,p,q,fn,fc);
		for(int time = 0; time < tmax; time++){
			if(time % tau == 0){
				code->introduceErrors();
			}
			code->getStarSyndromes();
			if(!code->hasSyndrome()){
				time = time + (tau - time % tau-1);
				continue;
			}
			this->simulationStep(code,automata,L);
			if(code->hasLogicalError()){
				code->getStarSyndromes();
				if(code->hasLogicalErrorAfterMWM()){
					#pragma omp critical
					{
						decayTimes.push_back(time);
						//std::cout << "Decay time: " << time << std::endl;
						//sum += time;
						//std::cout << "Average so far: " << sum/globalIndex << std::endl;
						//globalIndex++;
						//break;
						time = tmax;
					}
				}
			}
		}
		this->reset(code, automata, L);
	}
	return decayTimes;
}

void Simulator::init() {
	this->init(this->code, this->automata, this->L, this->Q, this->U, this->b, this->k,
			this->p, this->q, this->fn, this->fc);
}

void Simulator::init(ToricCode* &code, Automaton*** &automata, int L, int Q, int U, int b, int k,
		double p, double q, int fn, int fc){
	code = new ToricCode(L, p, q);
	automata = new Automaton**[L];
	for (int i = 0; i < L; i++) {
		automata[i] = new Automaton*[L];
		for (int j = 0; j < L; j++) {
			automata[i][j] = new Automaton(Point(j, L - i - 1), L, Q, U, b, k,p,fn,fc);
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
			automata[i][j]->localUpdates();
		}
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->executeFlips();
		}
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->signalUpdates();
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
			automata[i][j]->localUpdates();
		}
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->executeFlips();
		}
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			automata[i][j]->signalUpdates();
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
