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



void decayRun(const char* argv[]){
	int N,L,Q,k,U,b,tau,fn,fc;
	float p;
	bool printMovie = false;

	L=atoi(argv[1]);
	Q=atoi(argv[2]);
	k=atoi(argv[3]);
	U=atoi(argv[4]);
	N=atoi(argv[5]);
	tau=atoi(argv[6]);
	fn = atoi(argv[7]);
	fc = atoi(argv[8]);
	p = atof(argv[9]);
	b = isqrt(U);
	
	if(!printMovie){
		std::cout << "Running L = " << argv[1]
			<< ", Q = " << Q
			<< ", k = " << k
			<< ", U = " << U
			<< ", N = " << N
			<< ", tau = " << tau
			<< ", fn = " << fn
			<< ", fc = " << fc
			<< ", p = " << p
			<< std::endl;
	}


	//User defined failure probability
	std::vector<float> ps({p});
	
	
	//------- q = p----------
	//used for L=27
	/*std::vector<float> ps({0.00197,0.00184,0.00172});
		0.00161,0.0015,0.0014,0.00131,0.00123,
		0.00114,0.00107,0.001,0.00093,0.00087,
		0.00082,0.00076,0.00071});*/
	
	//used for L=3 and L=9
	/*std::vector<float> ps({0.0107,0.01,0.0093,0.0087,0.0082,
		0.0076,0.0072,0.0067,0.0062,0.0058,
		0.0054,0.0051,0.0047,0.0044,0.0041,0.0039,0.0036,
		0.0034,0.0032,0.00295,0.00276,0.00258,0.00241,
		0.00225,0.0021,0.00197,0.00184,0.00172,
		0.00161,0.0015,0.0014,0.00131,0.00123,
		0.00114,0.00107,0.001,0.00093,0.00087,
		0.00082,0.00076,0.00071});*/
	

	
	
	//------- q = 0----------

	// L=3 and L=9
	/*std::vector<float> ps({0.0076,0.0072,0.0067,0.0062,0.0058,
		0.0054,0.0051,0.0047,0.0044,0.0041,0.0039,0.0036,
		0.0034,0.0032,0.00295,0.00276,0.00258,0.00241,
		0.00225,0.0021,0.00197,0.00184,0.00172,
		0.00161,0.0015,0.0014});*?
	
	//--------rough test--------
	// L=3 and L=9
	/*std::vector<float> ps({0.0107,0.01,0.0093,0.0087,0.0082,
		0.0076,0.0072,0.0067,0.0062,0.0058,
		0.0054,0.0051,0.0047,0.0044,0.0041,0.0039,0.0036,
		0.0034,0.0032,0.00295,0.00276,0.00258,0.00241,
		0.00225,0.0021,0.00197,0.00184,0.00172,
		0.00161,0.0015,0.0014});*/
	// L=27
	/*std::vector<float> ps({0.00197,0.00184,0.00172,
		0.00161,0.0015,0.0014});*/
	
	//---------tau = infinite, q=0
	//used for L=3 and L=9
	/*std::vector<float> ps({0.051,0.047,0.044,0.041,0.039,0.036,
		0.034,0.032,0.0295,0.0276,0.0258,0.0241,
		0.0225,0.021,0.0197,0.0184,0.0172,
		0.0161,0.015});*/
	
	//---------tau = analysis, q=0 & q=p
	//L=9
	/*std::vector<float> ps({0.0107,0.01,0.0093,0.0087,0.0082,0.0076,0.0072,0.0067,0.0062,0.0058,
		0.0054,0.0051,0.0047,0.0044,0.0041,0.0039,0.0036,
		0.0034,0.0032,0.00295,0.00276,0.00258,0.00241,
		0.00225,0.0021,0.00197,0.00184,0.00172,
		0.00161,0.0015,0.0014,0.00131,0.00123,
		0.00114,0.00107,0.001,0.00093,0.00087,
		0.00082,0.00076,0.00071});*/
	
	
	
	// reverse the probability vector ps to do fast calculations for low p first
	std::reverse(ps.begin(),ps.end());
	int size = ps.size();
	int tmax = 1e9;

	Simulator sim;
	sim.setN(N);
	sim.setTmax(tmax);

	for(int i = 0; i<size; i++){
		std::list<int> data = sim.generateDecayCurve(L,Q,k,U,b,tau,ps[i],fn,fc,printMovie);
		if(!printMovie)
			cout << ps[i] << ",";
		int j = 0;
		int sum = 0;
		for(std::list<int>::const_iterator it = data.begin(); it != data.end(); it++){
			sum += *it;
		}
		if(!printMovie)
		cout << ((double)sum)/data.size() << ","
			<< average_error(data) << endl;
	}
}

void fnfcRun(const char* argv[]){
	int N,L,Q,k,U,b,tau;
	bool printMovie = false;

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


	int tmax = 1e9;
	list<int> frange;
	for(int i = 0; i <= U; i++){
		frange.push_back(i);
	}

	Simulator sim;
	sim.setN(N);
	sim.setTmax(tmax);

	for(double fn : frange){
		for(double fc : frange){
			list<int> data = sim.generateDecayCurve(L,Q,k,U,b,tau,p,fn,fc, printMovie);
			cout << fn << "," << fc << ",";
			int j = 0;
			int sum = 0;
			for(list<int>::const_iterator it = data.begin(); it != data.end(); it++){
				sum += *it;
			}
			cout << ((double)sum)/data.size() << ","
				<< average_error(data) << endl;
		}
	}
}

int main(int argc, const char* argv[]){
	switch(argc){
	case 6:
		// parameters: L, Q, k, U, N
		fnfcRun(argv);
		break;
	case 10:
		// parameters: L, Q, k, U, N, tau, fn, fc,  p
		decayRun(argv);
		break;
	}
	return 0;
}
