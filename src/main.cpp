/*
 * main.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <list>
#include "ToricCode.h"
#include "Location.h"
#include "Automaton.h"
#include "Simulator.h"
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>
#include <omp.h>

using namespace std;

/* function to calculate integer squareroot */
unsigned long isqrt(unsigned long x)
{
    register unsigned long op, res, one;

    op = x;
    res = 0;

    /* "one" starts at the highest power of four <= than the argument. */
    one = 1 << 30;  /* second-to-top bit set */
    while (one > op) one >>= 2;

    while (one != 0) {
        if (op >= res + one) {
            op -= res + one;
            res += one << 1;  // <-- faster than 2 * one
        }
        res >>= 1;
        one >>= 2;
    }
    return res;
}

/* function to calculate the error on the average of a given dataset */
double average_error(std::list<int> data) {
	int n = data.size();
    double mean=0.0, sum_deviation=0.0;
	int i;
	for(list<int>::const_iterator it = data.begin(); it != data.end(); it++)
		mean+=*it;
	mean=mean/n;
	for(list<int>::const_iterator it = data.begin(); it != data.end(); it++)
		sum_deviation+=(*it-mean)*(*it-mean);
	return sqrt(sum_deviation/n/(n-1));           
}

void thresholdRun(const char* argv[]){
	int N,L,Q,k,U,b,fn,fc, resolution;
	double pmin, pmax;

	pmin = 0.02;
	pmax = 0.05;
	resolution = 20;
	L=atoi(argv[1]);
	Q=atoi(argv[2]);
	k=atoi(argv[3]);
	U=atoi(argv[4]);
	N=atoi(argv[5]);
	fn = atoi(argv[6]);
	fc = atoi(argv[7]);

		b = isqrt(U);

	std::cout << "Running L = " << argv[1]
			<< ", Q = " << Q
			<< ", k = " << k
			<< ", U = " << U
			<< ", N = " << N
			<< ", fn = " << fn
			<< ", fc = " << fc
			<< std::endl;

	int tmax = 1e7;
	Simulator sim;
	sim.setN(N);
	sim.setTmax(tmax);
	sim.generateThresholdCurve(L,Q,k,U,b,resolution,pmin,pmax,fn,fc);
}

void decayRun(const char* argv[]){
	int N,L,Q,k,U,b,tau,fn,fc;

	L=atoi(argv[1]);
	Q=atoi(argv[2]);
	k=atoi(argv[3]);
	U=atoi(argv[4]);
	N=atoi(argv[5]);
	tau=atoi(argv[6]);
	fn = atoi(argv[7]);
	fc = atoi(argv[8]);
	b = isqrt(U);

	std::cout << "Running L = " << argv[1]
			<< ", Q = " << Q
			<< ", k = " << k
			<< ", U = " << U
			<< ", N = " << N
			<< ", tau = " << tau
			<< ", fn = " << fn
			<< ", fc = " << fc
			<< std::endl;


	std::vector<float> ps({0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01});
	//std::vector<float> ps({0.0007,0.0008});
	//std::vector<float> ps({0.0015,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01});
	//std::vector<float> ps({0.0015});
	
	// reverse the probability vector ps to do fast calculations for low p first
	std::reverse(ps.begin(),ps.end());
	int size = ps.size();
	int tmax = 1e7;

	Simulator sim;
	sim.setN(N);
	sim.setTmax(tmax);

	for(int i = 0; i<size; i++){
		std::list<int> data = sim.generateDecayCurve(L,Q,k,U,b,tau,ps[i],fn,fc);
		cout << ps[i] << ",";
		int j = 0;
		int sum = 0;
		for(std::list<int>::const_iterator it = data.begin(); it != data.end(); it++){
			sum += *it/tau;
		}
		cout << ((double)sum)/data.size() << ","
			<< average_error(data) << endl;
	}
}

void fnfcRun(const char* argv[]){
	int N,L,Q,k,U,b,tau;

	L=atoi(argv[1]);
	Q=atoi(argv[2]);
	k=atoi(argv[3]);
	U=atoi(argv[4]);
	N=atoi(argv[5]);
	b = isqrt(U);
	tau = 1;

	double p = 0.001;	
	std::cout << "Running L = " << L
			<< ", Q = " << Q
			<< ", k = " << k
			<< ", U = " << U
			<< ", N = " << N
			<< ", p = " << p
			<< std::endl;


	int tmax = 1e7;
	list<int> frange;
	for(int i = 0; i <= U; i++){
		frange.push_back(i);
	}

	Simulator sim;
	sim.setN(N);
	sim.setTmax(tmax);

	for(double fn : frange){
		for(double fc : frange){
			list<int> data = sim.generateDecayCurve(L,Q,k,U,b,tau,p,fn,fc);
			cout << fn << "," << fc << ",";
			int j = 0;
			int sum = 0;
			for(list<int>::const_iterator it = data.begin(); it != data.end(); it++){
				sum += *it/tau;
			}
			cout << ((double)sum)/data.size() << ","
				<< average_error(data) << endl;
		}
	}
}

void roughTestTestRun(const char* argv[]){
	long BG = 0;
	long BB = 0;
	long GB = 0;
	long GG = 0;
	long trials = atoi(argv[1]);
	float p = atof(argv[2]);
	#pragma omp parallel reduction(+:BG,BB,GB,GG)
	{
		int L = 5;
		ToricCode code(L,p,0.0);
		#pragma omp for
		for(auto i = 0; i < trials; i++){
			code.introduceErrors();
			code.getStarSyndromes();
			if(!code.hasLogicalError()){
				if(!code.hasLogicalErrorAfterMWM()){
					GG++;
				} else {
					GB++;
				}
			}
			else{
				if(!code.hasLogicalErrorAfterMWM()){
					BG++;
				} else {
					BB++;
				}
			}
			code.reset();
		}
	}
	cout << "GG: " << (double)GG/trials << endl;;
	cout << "GB: " << (double)GB/trials << endl;
	cout << "BG: " << (double)BG/trials << endl;
	cout << "BB: " << (double)BB/trials << endl;
}

int main(int argc, const char* argv[]){
	switch(argc){
	case 3:
		// parameters: trials, p
		roughTestTestRun(argv);
		break;
	case 6:
		// parameters: L, Q, k, U, N
		fnfcRun(argv);
		break;
	case 8:
		// parameters: L, Q, k, U, N, fn, fc
		thresholdRun(argv);
		break;
	case 9:
		// parameters: L, Q, k, U, N, tau, fn, fc
		decayRun(argv);
		break;
	}
	return 0;
}
