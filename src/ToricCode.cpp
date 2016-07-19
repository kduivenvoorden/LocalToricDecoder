/*
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "ToricCode.h"
#include "Qubit.h"
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <random>
#include "PerfectMatching.h"
#include <set>
#include <list>
#include <algorithm>

ToricCode::ToricCode(int L, double p, double q) : L(L), p(p), q(q), mt(generator()), distribution(0.0,1.0) {
	this->vertices = new int*[L];
	this->starSyndromes = new int*[L];
	this->faultySyndromes = new int*[L];
	for(int i = 0; i < L; i++){
		this->vertices[i] = new int[L];
		this->starSyndromes[i] = new int[L];
		this->faultySyndromes[i] = new int[L];
	}

	this->edges = new Qubit**[L];
	for(int i = 0; i < L; i++){
		this->edges[i] = new Qubit*[L];
		for(int j = 0; j < L; j++)
			this->edges[i][j] = new Qubit[2];
	}
	//initialize

	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			this->vertices[i][j] = 0;
			this->starSyndromes[i][j] = 0;
			this->faultySyndromes[i][j] = 0;
		}
	}
}

ToricCode::~ToricCode() {
	for(int i = 0; i < this->L; i++){
		delete[] this->vertices[i];
		delete[] this->starSyndromes[i];
		delete[] this->faultySyndromes[i];
	}
	delete[] this->vertices;
	delete[] this->starSyndromes;
	delete[] this->faultySyndromes;

	for(int i = 0; i < this->L; i++){
		for(int j = 0; j < this->L; j++){
			delete[] this->edges[i][j];
		}
		delete[] this->edges[i];
	}
	delete[] this->edges;
}

void ToricCode::setQ(double q){
	this->q = q;
}

/* the toriccode has periodic boundary conditions. therefore
 * vertex indices can in principle be arbitrary large. this
 * function is taking care of that case.
 */
int ToricCode::index(int i){
	if(i >= 0)
		return i % this->L;
	else
		return this->L-(abs(i) % this->L);
}

void ToricCode::reset(){
	for(int i = 0; i < this->L; i++){
		for(int j = 0; j < this->L; j++){
			this->vertices[i][j] = 0;
			this->starSyndromes[i][j] = 0;
			this->faultySyndromes[i][j] = 0;
			for(int k = 0; k < 2; k++){
				this->edges[i][j][k].bit = 0;
			}
		}
	}
}

/*
int ToricCode::starCheck(int i, int j){
	//double randnum = ((double) std::rand() / (RAND_MAX));
	//double randnum = this->distribution(this->mt);
	int ret = (this->edges[this->index(i-1)][this->index(j)][1].bit
			+ this->edges[this->index(i)][this->index(j-1)][0].bit
			+ this->edges[this->index(i)][this->index(j)][0].bit
			+ this->edges[this->index(i)][this->index(j)][1].bit)%2;
	return ret;
	//if (randnum < this->q)
		//return (ret + 1) % 2;
	//else
		//return ret;
}
*/

int ToricCode::starCheck(int i, int j){
	//double randnum = this->distribution(this->mt);
	int ret = (this->edges[this->index(i-1)][this->index(j)][1].bit
			+ this->edges[this->index(i)][this->index(j-1)][0].bit
			+ this->edges[this->index(i)][this->index(j)][0].bit
			+ this->edges[this->index(i)][this->index(j)][1].bit)%2
			^this->faultySyndromes[i][j];
	//int faulty = this->faultySyndromes[i][j];
	//if (randnum < this->q)
		//return (ret + 1) % 2;
	//else
	return ret;
}

int** ToricCode::getStarSyndromes(){
	for(int i = 0; i < this->L; i++){
		for(int j = 0; j < this->L; j++){
			this->starSyndromes[i][j] = this->starCheck(i,j);
		}
	}
	return starSyndromes;
}

void ToricCode::introduceErrors(){
	double randnum;
	for(int i = 0; i < this->L; i++){
		for(int j = 0; j < this->L; j++){
			for(int k = 0; k < 2; k++){
				//double randnum = ((double) std::rand() / (RAND_MAX));
				//randnum = ((double) std::rand() / (RAND_MAX));
				randnum = this->distribution(this->mt);
				if (randnum <= this->p)
					this->edges[i][j][k].flip();
			}
			randnum = this->distribution(this->mt);
			if (randnum <= this->q)
				this->faultySyndromes[i][j] = 1;
			else
				this->faultySyndromes[i][j] = 0;
		}
	}
}

bool ToricCode::hasSyndrome(){
	for(int i = 0; i < this->L; i++){
		for(int j = 0; j < this->L; j++){
			if(this->starSyndromes[i][j] == 1){
				return true;
			}
		}
	}
	return false;
}

/* function implementing the rough test */
bool ToricCode::hasLogicalError(Qubit*** edges, int L){
	std::vector<int> rowsums(L, 0);
	std::vector<int> columnsums(L, 0);
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			rowsums[i] += edges[i][j][1].bit;
		}
	}
	for(int j = 0; j < L; j++){
		for(int i = 0; i < L; i++){
			columnsums[j] += edges[i][j][0].bit;
		}
	}

	int rowsum = 0;
	int columnsum = 0;
	for(int i = 0; i < L; i++){
		rowsum += rowsums[i]%2;
		columnsum += columnsums[i]%2;
	}
	if(rowsum >= L/2 + 1 || columnsum >= L/2 + 1)
		return true;
	else
		return false;

}

bool ToricCode::hasLogicalError(){
	return this->hasLogicalError(this->edges, this->L);
}

void ToricCode::minimumWeightMatching(){
	int L = this->L;
	Qubit*** edges = this->edges;

	//setup copy of all edges to perform minimum weight matching
	std::vector<int> syndromeVerts;
	std::vector<std::vector<int> > syndromeVertsAddresses;
	std::vector<int> syndromeEdges;

	int numSyndromes = 0;
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			if(this->starSyndromes[i][j] == 1){
				syndromeVerts.push_back(numSyndromes);
				syndromeVertsAddresses.push_back({i,j});
				numSyndromes++;
			}
		}
	}

	int num_verts = numSyndromes;
	int num_edges = num_verts*(num_verts-1)/2;
	PerfectMatching *pm = new PerfectMatching(num_verts, num_edges);

	struct PerfectMatching::Options options;
	options.verbose = false;
	pm->options = options;

	for(int i = 0; i < num_verts; i++){
		for(int j = i+1; j < num_verts; j++){
			pm->AddEdge(syndromeVerts[i], syndromeVerts[j],
					this->getDistance(syndromeVertsAddresses[i][0], syndromeVertsAddresses[i][1],
							syndromeVertsAddresses[j][0], syndromeVertsAddresses[j][1]));
		}
	}
	pm->Solve();

	auto makeEdgeIndex = [this] (int i, int j){
		int max,min;
		if(i < j){
			max = j;
			min = i;
		} else {
			max = i;
			min = j;
		}
		return max*this->L*this->L + min;
	};
	std::list<int> doneFlips;
	for(int vert = 0; vert < num_verts; vert++){
		int partner = pm->GetMatch(syndromeVerts[vert]);
		int row1 = syndromeVertsAddresses[vert][0];
		int column1 = syndromeVertsAddresses[vert][1];
		int row2 = syndromeVertsAddresses[partner][0];
		int column2 = syndromeVertsAddresses[partner][1];
		int edgeIndex = makeEdgeIndex(vert,partner);
		if(std::find(doneFlips.begin(), doneFlips.end(), edgeIndex) != doneFlips.end()){
			continue;
		}
		//flip minimum weight path
		int horizontalFlipDirection = 0;
		int verticalFlipDirection = 0;

		int distance[4];
		int rowDistance = row1 - row2;
		int columnDistance = column1 - column2;
		int absRowDistance = abs(row1 - row2);
		int absColumnDistance = abs(column1 - column2);

		/*
		if(absRowDistance > this->L - absRowDistance){
			absRowDistance = this->L - absRowDistance;
			rowDistance = this->L - row1 - row2;
		}
		if(absColumnDistance > this->L - absColumnDistance){
			absColumnDistance = this->L - absColumnDistance;
			columnDistance = this->L - column1 - column2;
		}
		*/
		//int rowDistance[2];
		//int columnDistance[2];
		//rowDistance[0] = row1 - row2;
		//rowDistance[1] = row1 - (this->L - row2);
		//columnDistance[0] = column1 - column2;
		//columnDistance[1] = column1 - (this->L - column2);
		//int rowDistance = (row1 - row2);
		//int columnDistance = (column1 - column2);
		//distance[0] = absRowDistance + absColumnDistance;
		//distance[1] = this->L - absRowDistance + absColumnDistance;
		//distance[2] = absRowDistance + this->L - absColumnDistance;
		//distance[3] = this->L - rowDistance + this->L - absColumnDistance;
		distance[0] = absRowDistance + absColumnDistance;
		distance[1] = this->L - absRowDistance + absColumnDistance;
		distance[2] = absRowDistance + this->L - absColumnDistance;
		distance[3] = this->L - absRowDistance + this->L - absColumnDistance;
		int mindistance = *std::min_element(distance,distance+4);

		if(distance[0] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::S;
			else
				verticalFlipDirection = Location::N;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::E;
			else
				horizontalFlipDirection = Location::W;
		} else if(distance[1] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::N;
			else
				verticalFlipDirection = Location::S;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::E;
			else
				horizontalFlipDirection = Location::W;
		} else if(distance[2] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::S;
			else
				verticalFlipDirection = Location::N;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::W;
			else
				horizontalFlipDirection = Location::E;
		} else {
			if(rowDistance < 0)
				verticalFlipDirection = Location::N;
			else
				verticalFlipDirection = Location::S;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::W;
			else
				horizontalFlipDirection = Location::E;
		}

		if(absRowDistance > this->L - absRowDistance){
			absRowDistance = this->L - absRowDistance;
			rowDistance = this->L - row1 - row2;
		}
		if(absColumnDistance > this->L - absColumnDistance){
			absColumnDistance = this->L - absColumnDistance;
			columnDistance = this->L - column1 - column2;
		}

		Qubit* workingBit;
		for(int i = 1; i <= absRowDistance; i++){
			workingBit = this->getQubit(edges, row1, column1, verticalFlipDirection, i);
			workingBit->flip();
		}
		int offset;
		if(verticalFlipDirection == Location::N)
			offset = -absRowDistance;
		else
			offset = absRowDistance;
		for(int i = 1; i <= absColumnDistance; i++){
			workingBit = this->getQubit(edges, this->index(row1+offset), column1, horizontalFlipDirection, i);
			workingBit->flip();
		}
		doneFlips.push_back(edgeIndex);
	}
	delete pm;
}

bool ToricCode::hasLogicalErrorAfterMWM(){
	int L = this->L;
	//setup copy of all edges to perform minimum weight matching
	std::vector<int> syndromeVerts;
	std::vector<std::vector<int> > syndromeVertsAddresses;
	std::vector<int> syndromeEdges;

	int numSyndromes = 0;
	Qubit*** edges = new Qubit**[L];
	for(int i = 0; i < L; i++){
		edges[i] = new Qubit*[L];
		for(int j = 0; j < L; j++)
			edges[i][j] = new Qubit[2];
	}
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			for(int k = 0; k < 2; k++){
				edges[i][j][k].bit = this->edges[i][j][k].bit;
			}
			//syndromes here should not be faulty. correct for it
			if(this->starSyndromes[i][j]^faultySyndromes[i][j] == 1){
				syndromeVerts.push_back(numSyndromes);
				syndromeVertsAddresses.push_back({i,j});
				numSyndromes++;
			}
		}
	}

	int num_verts = numSyndromes;
	int num_edges = num_verts*(num_verts-1)/2;
	PerfectMatching *pm = new PerfectMatching(num_verts, num_edges);

	struct PerfectMatching::Options options;
	options.verbose = false;
	pm->options = options;

	for(int i = 0; i < num_verts; i++){
		for(int j = i+1; j < num_verts; j++){
			pm->AddEdge(syndromeVerts[i], syndromeVerts[j],
					this->getDistance(syndromeVertsAddresses[i][0], syndromeVertsAddresses[i][1],
							syndromeVertsAddresses[j][0], syndromeVertsAddresses[j][1]));
		}
	}
	pm->Solve();

	auto makeEdgeIndex = [this] (int i, int j){
		int max,min;
		if(i < j){
			max = j;
			min = i;
		} else {
			max = i;
			min = j;
		}
		//this->L squared added because there maybe maximum L squared vertices
		return max*this->L*this->L+ min;
	};
	std::list<int> doneFlips;
	for(int vert = 0; vert < num_verts; vert++){
		int partner = pm->GetMatch(syndromeVerts[vert]);
		int row1 = syndromeVertsAddresses[vert][0];
		int column1 = syndromeVertsAddresses[vert][1];
		int row2 = syndromeVertsAddresses[partner][0];
		int column2 = syndromeVertsAddresses[partner][1];
		int edgeIndex = makeEdgeIndex(vert,partner);
		if(std::find(doneFlips.begin(), doneFlips.end(), edgeIndex) != doneFlips.end()){
			continue;
		}
		//flip minimum weight path
		int horizontalFlipDirection = 0;
		int verticalFlipDirection = 0;

		int distance[4];
		int rowDistance = row1 - row2;
		int columnDistance = column1 - column2;
		int absRowDistance = abs(row1 - row2);
		int absColumnDistance = abs(column1 - column2);


		//int rowDistance[2];
		//int columnDistance[2];
		//rowDistance[0] = row1 - row2;
		//rowDistance[1] = row1 - (this->L - row2);
		//columnDistance[0] = column1 - column2;
		//columnDistance[1] = column1 - (this->L - column2);
		//int rowDistance = (row1 - row2);
		//int columnDistance = (column1 - column2);
		//distance[0] = absRowDistance + absColumnDistance;
		//distance[1] = this->L - absRowDistance + absColumnDistance;
		//distance[2] = absRowDistance + this->L - absColumnDistance;
		//distance[3] = this->L - rowDistance + this->L - absColumnDistance;
		distance[0] = absRowDistance + absColumnDistance;
		distance[1] = this->L - absRowDistance + absColumnDistance;
		distance[2] = absRowDistance + this->L - absColumnDistance;
		distance[3] = this->L - absRowDistance + this->L - absColumnDistance;
		int mindistance = *std::min_element(distance,distance+4);

		if(distance[0] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::S;
			else
				verticalFlipDirection = Location::N;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::E;
			else
				horizontalFlipDirection = Location::W;
		} else if(distance[1] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::N;
			else
				verticalFlipDirection = Location::S;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::E;
			else
				horizontalFlipDirection = Location::W;
		} else if(distance[2] == mindistance){
			if(rowDistance < 0)
				verticalFlipDirection = Location::S;
			else
				verticalFlipDirection = Location::N;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::W;
			else
				horizontalFlipDirection = Location::E;
		} else {
			if(rowDistance < 0)
				verticalFlipDirection = Location::N;
			else
				verticalFlipDirection = Location::S;
			if(columnDistance < 0)
				horizontalFlipDirection = Location::W;
			else
				horizontalFlipDirection = Location::E;
		}
		if(absRowDistance > this->L - absRowDistance){
			absRowDistance = this->L - absRowDistance;
			rowDistance = this->L - row1 - row2;
		}
		if(absColumnDistance > this->L - absColumnDistance){
			absColumnDistance = this->L - absColumnDistance;
			columnDistance = this->L - column1 - column2;
		}

		Qubit* workingBit;
		for(int i = 1; i <= absRowDistance; i++){
			workingBit = this->getQubit(edges, row1, column1, verticalFlipDirection, i);
			workingBit->flip();
		}
		int offset;
		if(verticalFlipDirection == Location::N)
			offset = -absRowDistance;
		else
			offset = absRowDistance;
		for(int i = 1; i <= absColumnDistance; i++){
			workingBit = this->getQubit(edges, this->index(row1+offset), column1, horizontalFlipDirection, i);
			workingBit->flip();
		}
		doneFlips.push_back(edgeIndex);
	}
	delete pm;
	//std::cout << "After minimum weight matching:" << std::endl;
	//this->dump(this->starSyndromes, edges, L);

	bool logicalError = this->hasLogicalError(edges, L);
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			delete[] edges[i][j];
		}
		delete[] edges[i];
	}
	delete[] edges;
	return logicalError;
}

int ToricCode::getDistance(int row1, int column1, int row2, int column2){
	int distance[4];
	int rowDistance = abs(row1 - row2);
	int columnDistance = abs(column1 - column2);
	distance[0] = rowDistance + columnDistance;
	distance[1] = this->L - rowDistance + columnDistance;
	distance[2] = rowDistance + this->L - columnDistance;
	distance[3] = this->L - rowDistance + this->L - columnDistance;
	return *std::min_element(distance,distance+4);
}

Qubit** ToricCode::getQubits(int i, int j){
	Qubit** ret = new Qubit*[4];
	ret[0] = &(this->edges[this->index(i-1)][this->index(j)][1]);
	ret[1] = &(this->edges[this->index(i)][this->index(j-1)][0]);
	ret[2] = &(this->edges[this->index(i)][this->index(j)][0]);
	ret[3] = &(this->edges[this->index(i)][this->index(j)][1]);
	return ret;
}

Qubit* ToricCode::getQubit(Qubit*** edges, int i, int j, int l, int multiplicity = 1){
	switch(l){
	case Location::N:
		return &(edges[this->index(i-multiplicity)][this->index(j)][1]);
	case Location::W:
		return &(edges[this->index(i)][this->index(j-multiplicity)][0]);
	case Location::E:
		return &(edges[this->index(i)][this->index(j+multiplicity-1)][0]);
	case Location::S:
		return &(edges[this->index(i+multiplicity-1)][this->index(j)][1]);
	default:
		//TODO: pass some kind of error message
		return edges[0][0];
	}
}

void ToricCode::flip(Qubit*** edges, int i, int j, int l){
	switch(l){
	case Location::E:
		edges[this->index(i)][this->index(j)][0].flip();
		break;
	case Location::S:
		edges[this->index(i)][this->index(j)][1].flip();
		break;
	case Location::N:
		edges[this->index(i-1)][this->index(j)][1].flip();
		break;
	case Location::W:
		edges[this->index(i)][this->index(j-1)][0].flip();
		break;
	}
}


void ToricCode::flip(int i, int j, int l){
	this->flip(this->edges, i, j, l);
}

void ToricCode::dump(int** starSyndromes, Qubit*** edges, int L){
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			if (this->starSyndromes[i][j] == 1)
				std::cout << "s";
			else
				std::cout << "*";
			if (edges[i][j][0].bit == 1)
				std::cout << "-";
			else
				std::cout << " ";
		}
		std::cout << std::endl;
		for(int j = 0; j < L; j++){
			if (edges[i][j][1].bit == 1)
				std::cout << "| ";
			else
				std::cout << "  ";
		}
		std::cout << std::endl;
	}
}
void ToricCode::dump(){
	this->dump(this->starSyndromes, this->edges, this->L);
}
